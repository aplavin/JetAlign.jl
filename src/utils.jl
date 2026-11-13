# A memoized derived value computed on a background task: recomputes only when `inputs` change,
# coalescing while busy (at most one computation in flight). The UI never blocks; it keeps showing
# the last published result. Inputs are compared with `==`.
mutable struct AsyncJob
    task::Union{Task,Nothing}
    inputs                                  # inputs of the latest started computation
    published                               # (; result, inputs) or nothing
end
AsyncJob() = AsyncJob(nothing, nothing, nothing)

# Call every frame with the current `inputs`; returns the latest `(; result, inputs)`, or nothing.
function tick!(compute, job::AsyncJob, inputs)
    if !isnothing(job.task) && istaskdone(job.task)
        try
            job.published = (; result=fetch(job.task), inputs=job.inputs)
        catch e
            @error "background computation failed" exception=(e, catch_backtrace())
        end
        job.task = nothing
    end
    if isnothing(job.task) && inputs != job.inputs
        job.inputs = inputs
        job.task = Threads.@spawn compute(inputs)
    end
    job.published
end

# Memoize `f()` on `deps`, compared with `===`: identity for arrays, value for isbits.
mutable struct CacheCell
    deps
    value
end
CacheCell() = CacheCell(nothing, nothing)

function cached(f, deps, c::CacheCell)
    if c.deps !== deps
        c.value = f()
        c.deps = deps
    end
    c.value
end


# UI state persisted per file pair: a Serialization blob of the state's properties in a SQLite table.
struct StateStore{D,K}
    dict::D
    key::K
end
open_store(dbfile, files) = let key = @NamedTuple{f1::String, f2::String}(NTuple{2}(files))
    StateStore(SQLDictionary{typeof(key), @NamedTuple{ui::Any}}(SQLite.DB(dbfile), "ui_state"), key)
end

"The stored state of `store`, or `default` if none is stored or its type differs from `default`'s."
restore(::Nothing, default) = default
function restore(store::StateStore, default)
    stored = get(store.dict, store.key, nothing)
    isnothing(stored) && return default
    props = getproperties(default)
    if typeof(stored.ui) != typeof(props)
        @warn "Previous UI state is stored, but types do not match. Resetting to default."
        return default
    end
    setproperties(default, stored.ui)
end

save!(::Nothing, state) = nothing
save!(store::StateStore, state) = set!(store.dict, store.key, (; ui=getproperties(state)))


"Shifts stored in `file` for the image pair `files`; errors if `file` was saved for other files."
function load_shifts(file, files)
    data = TOML.parsefile(file)
    data["files"] == files || error("Shifts file $file is for files $(data["files"]), not for $files")
    map(SVector{2,Float64}, data["shifts"])
end

save_shifts(file, files, shifts) = open(io -> TOML.print(io, Dict("files" => files, "shifts" => shifts)), file, "w")
