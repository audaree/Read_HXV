## Julia program to read a selected .HXV file and display 30-minute time series plots
## JW December 2022
#using ContinuousWavelets 
using CSV
using Dates, DataFrames, Distributions, DSP
using FFTW
using Gtk
using LaTeXStrings
using NativeFileDialog
using Plots
using Printf
using Statistics #, StatsPlots
using Tk


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
    
    
function calc_csd_welch(X,Y,sample_frequency,len,olap)
#####################################

    segmentsX = arraysplit(X,len,olap)
    segmentsY = arraysplit(Y,len,olap)

#    N = length(X)
    dt = 1/sample_frequency
    
    csd = []

    for i in eachindex(segmentsX)
        
        x = fftshift(fft(segmentsX[i] )) 
        y = fftshift(fft(segmentsY[i] ))

        push!(csd, (2 * dt^2 / len) .* (x .* conj(y)) .* sample_frequency .* tukey)

    end
        
    return(mean(csd,dims=1)[1])
    
end    # calc_csd_welch()

        
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


################################################
################################################
##           START OF MAIN PROGRAM
################################################
################################################

# Widen screen for better viewing
##display("text/html", "<style>.container { width:100% !important; }</style>")

hxv_directory = pick_folder()

# build list of all hxv files in selected directory
hxv_files = filter(x->occursin(".hxv",x), readdir(hxv_directory));
hxv_files = hxv_files[findall(x->endswith(uppercase(x), ".HXV"), hxv_files)];

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
    global wse_df = DataFrame(Date = unix2datetime.(datetime2unix.(start_date) .+ (0:1/1.28:1800-1/1.28)), Heave = zeros(2304), North = zeros(2304), West = zeros(2304));

    # populate the df based on sequence numbers
    wse_df = calc_wse(infil, wse_df, start_date)
    is_gps = process_spectrum_file()

    N = length(wse_df.Heave)
    sample_frequency = 1.28
    dt = 1/sample_frequency
    global tukey = DSP.Windows.tukey(256,66/256)

     csd_vv = calc_csd_welch(wse_df.Heave,wse_df.Heave,sample_frequency,256,0)
    csd_nn = calc_csd_welch(wse_df.North,wse_df.North,sample_frequency,256,0)
    csd_ww = calc_csd_welch(wse_df.West,wse_df.West,sample_frequency,256,0)
    csd_nw = calc_csd_welch(wse_df.North,wse_df.West,sample_frequency,256,0)
    csd_vn = calc_csd_welch(wse_df.Heave,wse_df.North,sample_frequency,256,0)
    csd_vw = calc_csd_welch(wse_df.Heave,wse_df.West,sample_frequency,256,0)

    f2 = fftshift(fftfreq(length(csd_vv), 1/dt))
    f2_pos_vals = findall( x -> x >= 0, f2 )
    
    # Calculate the co- and quad- spectra
    cvv = real(csd_vv[f2_pos_vals])
    cnn = real(csd_nn[f2_pos_vals])
    cww = real(csd_ww[f2_pos_vals])
    cnw = real(csd_nw[f2_pos_vals])
    qvn = imag(csd_vn[f2_pos_vals])
    qvw = imag(csd_vw[f2_pos_vals])

    # Calculate the Fourier Coefficients
    a1 = qvn ./ (cvv .* (cnn .+ cww)).^0.5
    b1 = qvw ./ (cvv .* (cnn .+ cww)).^0.5

    a2 = (cnn .- cww) ./ (cnn .+ cww)
    b2 = (2 .* cnw) ./ (cnn .+ cww);

    # check for any Nans in the data
    replace!(a1, NaN=>0)
    replace!(b1, NaN=>0)
    replace!(a2, NaN=>0)
    replace!(b2, NaN=>0);

    # Calculate the Centred Fourier Coefficients
    θ₀ = atan.(b1,a1)
    m1 = (a1.^2 .+ b1.^2).^0.5
    m2 = a2.*cos.(2*θ₀) .+ b2.*sin.(2*θ₀)
    n2 = -a2.*sin.(2*θ₀) .+ b2.*cos.(2*θ₀)

    # Calculate the spread
    σc = (2 .* (1 .- m1)).^0.5;

    N = length(wse_df.Heave)
    f1 = 0.005:0.005:0.64

    p1 = plot(f1,a1, lw=2, c=:blue, label="a₁", fg_legend = :false, title=("a₁ and b₁"))
    p1 = plot!(f1,b1, lw=2, c=:red, label="b₁")

    p2 = plot(f1,a2, lw=2, c=:blue, label="a₂", fg_legend = :false, title=("a₂ and b₂"))
    p2 = plot!(f1,b2, lw=2, c=:red, label="b₂")

    p3 = plot(f1,m1, lw=2, c=:blue, label="m₁", fg_legend = :false, title=("m₁ and m₂"))
    p3 = plot!(f1,m2, lw=2, c=:red, label="m₂")

    p4 = plot(f1,n2, lw=2, c=:blue, label="n₂", fg_legend = :false, title=("n₂"))

    coefficients_plot = plot(p1, p2, p3, p4, layout = (2,2), framestyle = :box, titlefontsize=10, xlabel="Frequency (Hz)", xlabelfontsize=9,
            leftmargin = 12Plots.mm, bottommargin = 12Plots.mm, grid=true, size=(1400, 800), colorbar=false, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1)

    display(coefficients_plot)

    direction = (rad2deg.(π/2 .- atan.(b1,a1)) .+ 450).%360
    dir1 = direction
    if is_gps
        dir_type = " True"
    else
        dir_type = " Magnetic"
    end
    
    p1 = plot(f2[f2_pos_vals], direction, lw=3, c=:blue, yflip=true, ylim=[0,360], yticks = 0:30:360, grid=true, label="Direction", xlabel="Frequency (Hz)",
        ylabel="Direction ("*L"^o"*dir_type*")", legend=:topright)
    p1 = plot!(f2[f2_pos_vals], direction .+ rad2deg.(σc), fillrange = direction .- rad2deg.(σc), fillalpha = 0.25, lw=0, c = 1, label = "Spread")
    
    max_spec = max(maximum(power(welch_pgram(wse_df.Heave, 256, 0; fs=sample_frequency, window=tukey))),maximum(real(csd_vv[f2_pos_vals]))) * 1.05
    
    p1 = plot!(twinx(),freq(welch_pgram(wse_df.Heave, 256, 0; fs=sample_frequency, window=tukey)),ylabel="Spectral Density (m"*L"^2"*"/Hz.)",
        power(welch_pgram(wse_df.Heave, 256, 0; fs=sample_frequency, window=tukey)),c=:red,lw=3, label="DSP Welch\n", ylim=[0,max_spec], legend=:bottomright)
    p1 = plot!(twinx(),f2[f2_pos_vals], real(csd_vv[f2_pos_vals]), lw=2, c=:yellow, label="Welch calc. JW", ylim=[0,max_spec], legend=:bottomright)
    
    p_all = plot(p1, framestyle = :box, xlim=[0,0.64], title="Spectra and Direction",
            leftmargin = 15Plots.mm, rightmargin = 25Plots.mm, bottommargin = 15Plots.mm, size=(1400, 600), colorbar=false, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1, fg_legend = :false)
    
    display(p_all)
end