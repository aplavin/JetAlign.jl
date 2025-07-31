using TestItems
using TestItemRunner
@run_package_tests


@testitem "imalign_interactive smoke test" begin
    using JetAlign
    using JetAlign.VLBIFiles
    using JetAlign.Makie
    using JetAlign.AccessorsExtra
    
    # Get test data from VLBIFiles package
    vlbifiles_test_dir = joinpath(pkgdir(VLBIFiles), "test", "data")
    map_fits = joinpath(vlbifiles_test_dir, "map.fits")
    
    # Test that imalign_interactive runs without error (don't display the figure)
    imalign_interactive([map_fits, map_fits]; display=false)
    Makie.colorbuffer(current_figure())
    imalign_interactive([map_fits, map_fits]; display=false)
    Makie.colorbuffer(current_figure())

    @showtime imalign_interactive([map_fits, map_fits]; display=false, persiststate=false)
    fig = current_figure()
    @showtime Makie.colorbuffer(fig)

    sliders = @getall fig |> RecursiveOfType(Makie.Slider)
    checkboxes = @getall fig |> RecursiveOfType(Makie.Checkbox)
    @test length(sliders) > 0
    # @test length(checkboxes) > 0
    for _ in 1:30
        for cb in checkboxes
            cb.checked[] = rand(Bool)
        end
        for sl in sliders
            sl.selected_index[] = rand(eachindex(sl.range[]))
        end
        Makie.colorbuffer(fig)
    end
end

@testitem "align_strings" begin
    using JetAlign: align_strings
    
    # Test exact match
    @test align_strings("hello", "hello") == [(:match, 'h', 'h'), (:match, 'e', 'e'), (:match, 'l', 'l'), (:match, 'l', 'l'), (:match, 'o', 'o')]
    
    # Test substitution
    @test align_strings("hello", "hallo") == [(:match, 'h', 'h'), (:substitute, 'e', 'a'), (:match, 'l', 'l'), (:match, 'l', 'l'), (:match, 'o', 'o')]
    
    # Test insertion
    @test align_strings("cat", "cart") == [(:match, 'c', 'c'), (:match, 'a', 'a'), (:insert, nothing, 'r'), (:match, 't', 't')]
    
    # Test deletion
    @test align_strings("testing", "test") == [(:match, 't', 't'), (:match, 'e', 'e'), (:match, 's', 's'), (:match, 't', 't'), (:delete, 'i', nothing), (:delete, 'n', nothing), (:delete, 'g', nothing)]
    
    # Test empty strings
    @test align_strings("", "") == []
    
    # Test one empty string
    @test align_strings("abc", "") == [(:delete, 'a', nothing), (:delete, 'b', nothing), (:delete, 'c', nothing)]
    
    @test align_strings("", "xyz") == [(:insert, nothing, 'x'), (:insert, nothing, 'y'), (:insert, nothing, 'z')]
    
    # Test complex case with multiple operations
    @test align_strings("kitten", "sitting") == [(:substitute, 'k', 's'), (:match, 'i', 'i'), (:match, 't', 't'), (:match, 't', 't'), (:substitute, 'e', 'i'), (:match, 'n', 'n'), (:insert, nothing, 'g')]
    
    # Test Cyrillic characters
    @test align_strings("привет", "привет") == [(:match, 'п', 'п'), (:match, 'р', 'р'), (:match, 'и', 'и'), (:match, 'в', 'в'), (:match, 'е', 'е'), (:match, 'т', 'т')]
    
    @test align_strings("кот", "кит") == [(:match, 'к', 'к'), (:substitute, 'о', 'и'), (:match, 'т', 'т')]
    
    @test align_strings("мир", "мира") == [(:match, 'м', 'м'), (:match, 'и', 'и'), (:match, 'р', 'р'), (:insert, nothing, 'а')]
    
    # Test mixed scripts
    @test align_strings("hello мир", "hello мира") == [(:match, 'h', 'h'), (:match, 'e', 'e'), (:match, 'l', 'l'), (:match, 'l', 'l'), (:match, 'o', 'o'), (:match, ' ', ' '), (:match, 'м', 'м'), (:match, 'и', 'и'), (:match, 'р', 'р'), (:insert, nothing, 'а')]
    
    # Test other Unicode characters (emoji, mathematical symbols)
    @test align_strings("😀😃", "😀😄") == [(:match, '😀', '😀'), (:substitute, '😃', '😄')]
    
    @test align_strings("α + β", "α - β") == [(:match, 'α', 'α'), (:match, ' ', ' '), (:substitute, '+', '-'), (:match, ' ', ' '), (:match, 'β', 'β')]
    
    # Test Chinese/Japanese characters
    @test align_strings("你好", "你好") == [(:match, '你', '你'), (:match, '好', '好')]
    
    @test align_strings("こんにちは", "こんばんは") == [(:match, 'こ', 'こ'), (:match, 'ん', 'ん'), (:substitute, 'に', 'ば'), (:substitute, 'ち', 'ん'), (:match, 'は', 'は')]
end

@testitem "_" begin
    import Aqua
    Aqua.test_all(JetAlign; piracies=(;broken=true))

    import CompatHelperLocal as CHL
    CHL.@check()
end
