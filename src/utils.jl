round_(x::Number, precision::Number) = round(Integer, x / precision) * precision

function setup_state_persistence(statekey, stateobs_single)
    state_db_file = joinpath(@get_scratch!("cache"), "cache.db")
    state_db = SQLite.DB(state_db_file)
    state_dict = SQLDictionary{typeof(statekey), @NamedTuple{ui::Any}}(state_db, "ui_state")
    stored_state = get(state_dict, statekey, nothing)
    if !isnothing(stored_state)
        if typeof(stored_state.ui) == typeof(stateobs_single[])
            stateobs_single[] = stored_state.ui
        else
            @warn "Previous UI state is stored, but types do not match. Resetting to default."
        end
    end
    on(Makie.Observables.throttle(2, stateobs_single)) do _
        set!(state_dict, statekey, (;ui=stateobs_single[]))
    end
end
