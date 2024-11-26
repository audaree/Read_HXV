## Julia program to locate GPS errors in .HXV files
## JW January 2003

#using ContinuousWavelets 
using CSV
using Dates, DataFrames, Distributions, DSP
using Gtk
using LaTeXStrings
using NativeFileDialog
using Plots
using Printf
##using Statistics #, StatsPlots
using Tk
##import Pkg; Pkg.add("CurveFit")
using  CurveFit

function get_displacements(arry)
#####################################
    
    displacements = []

    if length(arry[1]) == 3
    
        for i in arry
            append!(displacements,parse(Int, SubString.(i, 1, 1), base=16)*16^2 + parse(Int, SubString.(i, 2, 2), base=16)*16^1 + parse(Int, SubString.(i, 3, 3), base=16)*16^0)
        end
        
    else
        
        for i in arry
            append!(displacements,parse(Int, SubString.(i, 1, 1), base=16)*16^1 + parse(Int, SubString.(i, 2, 2), base=16)*16^0)
        end
        
    end

    displacements[findall(>=(2048), displacements)] = 2048 .- displacements[findall(>=(2048), displacements)];
    
    return(displacements./100)
    
    end     # get_displacements()


function get_HNW(infil)
#####################################
        
    global df = DataFrame(CSV.File(infil,header=0, delim=",", types=String));

    # Calculate sequence numbers
    arry = SubString.(df.Column1, 3, 4)

    global sequence = []

    for i in arry
        append!(sequence,parse(Int, SubString.(i, 1, 1), base=16)*16^1 + parse(Int, SubString.(i, 2, 2), base=16)*16^0)
    end

    # Calculate heave WSEs
    arry = SubString.(df.Column3, 1, 3);
    heave = get_displacements(arry);

    # Calculate north WSEs
    arry = SubString.(df.Column3, 4, ) .* SubString.(df.Column4, 1, 2)
    north = get_displacements(arry);

    # Calculate north WSEs
    arry = SubString.(df.Column4, 3, 4) .* SubString.(df.Column5, 1, 1)
    west = get_displacements(arry);

    return(heave, north, west)

    end    # get_HNW()


function calc_wse(infil, wse_df, start_date)
#####################################    
    
    heave, north, west = get_HNW(infil)
    
    # Identify any gaps in the recorded data
    tt = [0]
    append!(tt,diff(sequence))
    tt[tt.<0] .+= 256;
    tt1 = cumsum(tt);
    
    if length(tt1) > 2304
        tt1 = tt1[1:2304]
    end

    [wse_df[tt1[i]+1,2] = heave[i] for i in eachindex(tt1)];
    [wse_df[tt1[i]+1,3] = north[i] for i in eachindex(tt1)];
    [wse_df[tt1[i]+1,4] = west[i] for i in eachindex(tt1)];
    
    return(wse_df)
    
    end    # calc_wse()


function process_spectrum_file()
    global is_gps = false
    sync_word_location = findall(x -> x == "7FFF", df.Column2)

    for j in sync_word_location

        i = df.Column2[j+1]
        word_number = parse(Int, SubString.(i, 1, 1), base=16)*16^0
        word = parse(Int, SubString.(i, 2, 2), base=16)*16^2 + parse(Int, SubString.(i, 3, 3), base=16)*16^1 + parse(Int, SubString.(i, 4, 4), base=16)*16^0
    #    println(j,' ',word_number,' ',word)

        # Test whether buoy is MkIII or DWR-G - see p.51 Table 5.7.5a. Organization and significance of the system file data 
        # If DWR-G:
        #     Av0 = 0; Ax0 = 0; Ay0 = 0; O = 0; and Inclination = 0
        
        if (word_number == 7 && word == 0)
            is_gps = true
        end

    end
    
    return(is_gps)
    
    end    # process_spectrum_file()


function show_gps_error(df, wse_df, datawell_filter_df)
################################################
    
    # get a list of the gps errors in the selected record
    gps_errors = findall(isodd,parse.(Int,SubString.(string.(df.Column4), 2, 2), base = 16))

    if isempty(gps_errors)
        
        println("No GPS errors in this record")
        
    else
        
        println(length(gps_errors)," GPS errors in this record")
    
        # start at the first gps error
        i=1

        #Identify where GPS error is centered in wse
        center = gps_errors[i]


        # add df column of filter points from start of record (values less than zero will not be included in error correction process)
        datawell_filter_df.Location = datawell_filter_df.Points .+ center

        # identify whether parts of the filter extend beyond start or end of record - if so, they will be truncated
        visible_filter_points = findall(x-> start_point < x < end_point, datawell_filter_df.Points .+ center)

        # 30-minute record of WSE's
        p1 = Plots.plot(wse_df.Date, wse_df.Heave, c=:lightgrey, label="Heave") 

        # Section of WSE's centred on GPS error spanning length of Datawell filter
        p1 = Plots.plot!(wse_df.Date[datawell_filter_df.Location[visible_filter_points]], wse_df.Heave[datawell_filter_df.Location[visible_filter_points]], 
            c=:blue, alpha=0.5, label="Heave with GPS error") 

        # Datawell filter
        p1 = Plots.plot!(wse_df.Date[datawell_filter_df.Location[visible_filter_points]],datawell_filter_df.Column1[visible_filter_points], 
            c=:red, label="Datawell GPS filter") 

        # Location of GPS error flag
        p1 = vline!([wse_df.Date[center]], c=:yellow, lw=2, ls=:dash,label="GPS error location")

        title_string = Dates.format(wse_df.Date[1], "dd-mm-yyyy HH:MM")
        
        # display plots to screen
        tm_tick = range(first(wse_df.Date),last(wse_df.Date),step=Minute(5))
        ticks = Dates.format.(tm_tick,"MM:SS")

        show_gps_errors = Plots.plot(p1, size = (1600, 600), title=title_string, titlefontsize=10, framestyle = :box, fg_legend=:transparent, legend=:bottomleft,
            grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1,xlim=(first(wse_df.Date),last(wse_df.Date)), xticks=(tm_tick,ticks))

        savefig(".\\output_plot.png")

        display(show_gps_errors)
        
    end
    
    return()
    
    end    # show_gps_error()


function fix_gps_errors(wse_df,datawell_filter_df)
# function to apply polynomial fit to WSE's affected by GPS errors
# uses selectable offset value to fine-tune result
##########################################################    
    
    # make a 'deep' copy of the heave data
    heave = deepcopy(wse_df.Heave)

    # locate GPS errors
    gps_errors = findall(isodd,parse.(Int,SubString.(string.(df.Column4), 2, 2), base = 16))
    
    if isempty(gps_errors)
        
        println("No GPS errors in this record")
        
    else
        
        println(length(gps_errors)," GPS errors in this record")
        
        for ii in gps_errors

            error_center = ii

            # add df column of filter points from start of record (values less than zero will not be included in error correction process)
            datawell_filter_df.Location = datawell_filter_df.Points .+ error_center;

            # User-selected offset either side of GPS error
            lower_offset = upper_offset = 50

            if error_center <= lower_offset
                lower_offset = 2
            end

            if error_center+upper_offset > 2304
                upper_offset = 2304 - error_center
            end

            aa = error_center-lower_offset:error_center+upper_offset
            bb = findall(x-> -lower_offset <= x <= upper_offset, datawell_filter_df.Points)

            # Fit curve to subset of heave before GPS error
            cc = error_center-lower_offset:error_center
            fit1 = curve_fit(Polynomial, cc, heave[cc], 2)
            yfit1 = fit1.(cc)
            yfit1[length(yfit1)] = 0.0

            # Fit curve to subset of heave after GPS error
            dd = error_center:error_center+upper_offset
            fit2 = curve_fit(Polynomial, dd, heave[dd], 2)
            yfit2 = fit2.(dd)
            yfit2[1] = 0.0

            # apply polynomial results to wse's on both sides of GPS error
            heave[cc] .= heave[cc]-yfit1
            heave[dd] .= heave[dd]-yfit2
            heave[ii] = 0.0    # set wse at GPS error location to 0

        end
    
    end
    
    return(heave)
    
    end     # fix_gps_errors()


################################################
################################################
##           START OF MAIN PROGRAM
################################################
################################################

# Widen screen for better viewing
display("text/html", "<style>.container { width:100% !important; }</style>")

hxv_directory = pick_folder()

###############################################################################
# build list of all hxv files in selected directory
hxv_files = filter(x->occursin(".hxv",x), readdir(hxv_directory));
hxv_files = hxv_files[findall(x->endswith(uppercase(x), ".HXV"), hxv_files)];

# read Datawell filter data from .csv file to df
GPS_file = ".\\Datawell_GPS_filter.csv"
datawell_filter_df = CSV.read(GPS_file, DataFrame, header=false)

# add df column of filter point numbers (should be -378 to +378, centered on 0)
filter_range = Int(trunc(nrow(datawell_filter_df)/2))
filter_points = -filter_range:filter_range
insertcols!(datawell_filter_df,1,:Points =>filter_points)
###############################################################################

# need to specify range over which filter can extend
start_point = 0
end_point = 2304

w = Toplevel("Select Date", 235, 600)
tcl("pack", "propagate", w, false)
f = Frame(w)
pack(f, expand=true, fill="both")

f1 = Frame(f)
lb = Treeview(f1, hxv_files)
scrollbars_add(f1, lb)
pack(f1,  expand=true, fill="both")

tcl("ttk::style", "configure", "TButton", foreground="blue", font="arial 16 bold")
b = Button(f, "Ok")
pack(b)

bind(b, "command") do path
    
    global file_choice = get_value(lb);   

    # Select a HXV file
    global infil = hxv_directory * "\\" * file_choice[1]
    println("Selected ",infil)

    # extract the datetime from the file name
    date_str = split(infil,".")[1]
    ll = length(date_str)
    start_date = DateTime.(date_str[ll-16:ll-1], "yyyy-mm-ddTHHhMMZ")

    # create df of 2304 rows, each 0.78s apart
    global wse_df = DataFrame(Date = unix2datetime.(datetime2unix.(start_date) .+ (0:1/1.28:1800-1/1.28)),
        Heave = zeros(2304), North = zeros(2304), West = zeros(2304));

    # populate the df based on sequence numbers
    wse_df = calc_wse(infil, wse_df, start_date)
    is_gps = process_spectrum_file()
    
##    show_gps_error(df, wse_df,datawell_filter_df)
    fixed_heave = fix_gps_errors(wse_df,datawell_filter_df)
    
    # plot of initial heave with GPS errors
    p1 = Plots.plot(wse_df.Date,wse_df.Heave, c=:yellow, lw=2,  label="Uncorrected heave")

    # plot of correction applied to initial heave
    #p1 = plot!(wse_df.Date,wse_df.Heave-fixed_heave,lw=2)

    gps_errors = findall(isodd,parse.(Int,SubString.(string.(df.Column4), 2, 2), base = 16))

    # plot of corrected Heave
    p1 = Plots.plot!(wse_df.Date,fixed_heave, c=:blue, label="Filtered Heave")

    # plot of gps errors
    """
    for i in gps_errors
        p1 = vline!([wse_df.Date[i]], lw=0.5, ls=:dash, c=:red, label="")
    end
    """

    title_string = Dates.format(wse_df.Date[1], "dd-mm-yyyy HH:MM")

    # display plots to screen
    tm_tick = range(first(wse_df.Date), last(wse_df.Date), step=Minute(1))
    ticks = Dates.format.(tm_tick,"MM:SS")

    fix_gps = Plots.plot(p1,size = (1600, 600), title=title_string, titlefontsize=10, framestyle = :box, fg_legend=:transparent, legend=:bottomleft,
                    xlim=(first(wse_df.Date),last(wse_df.Date)), grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1,xticks=(tm_tick,ticks))

    savefig(".\\output_plot.png")

    display(fix_gps)
    
end