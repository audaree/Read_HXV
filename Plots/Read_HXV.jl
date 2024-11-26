## Julia program to read a selected .HXV file and display 30-minute time series plots
## JW December 2022
#using ContinuousWavelets 
using CSV
using Dates, DataFrames, Distributions, DSP
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


function spike_value(wse)
#####################################
    
    median_value = median(wse)
    std_value = std(wse)

    return(median_value + 3*std_value)

    end    # spike_value()


function plot_wses(wse_df, is_gps)
#####################################
    
    spike = spike_value(wse_df.Heave)
    heave_spikes = findall(i->(i>=spike), abs.(wse_df.Heave));

    spike = spike_value(wse_df.North)
    north_spikes = findall(i->(i>=spike), abs.(wse_df.North));

    spike = spike_value(wse_df.West)
    west_spikes = findall(i->(i>=spike), abs.(wse_df.West));
    
    # create plots of heave, north, and west
    title_string = Dates.format(first(wse_df.Date), "dd/mm/yyyy HH:MM") # * " UTC"
    p1_hnw = scatter(wse_df[heave_spikes,:].Date, wse_df[heave_spikes,:].Heave, label="", ylabel="Heave", markershape=:circle, ms=4, mc=:white, ma=1, msc=:red, msa=0.25, msw=0.5)
    p1_hnw = plot!(wse_df.Date,wse_df.Heave, label="", c="#4a536b", lw=0.5, title=title_string, titlefontsize=12) ##last(split(infil,"\\")))

    # get plotting limits
    x_lim1 = xlims(p1_hnw)[1]; y_lim1 = ylims(p1_hnw)[1]
    x_lim2 = xlims(p1_hnw)[2]; y_lim2 = ylims(p1_hnw)[2]

    p2_hnw = scatter(wse_df[north_spikes,:].Date, wse_df[north_spikes,:].North, label="", ylabel="North", markershape=:circle, ms=4, mc=:white, ma=1, msc=:red, msa=0.25, msw=0.5)
    p2_hnw = plot!(wse_df.Date,wse_df.North, label="", c="#aed6dc", lw=0.5)
    p3_hnw = scatter(wse_df[west_spikes,:].Date,wse_df[west_spikes,:].West, label="", ylabel="West", markershape=:circle, ms=4, mc=:white, ma=1, msc=:red, msa=0.25, msw=0.5)
    p3_hnw = plot!(wse_df.Date,wse_df.West, label="", c="#ff9a8d", lw=0.5)

    hline!(p1_hnw, [0], lw=0.5, label="")
    hline!(p2_hnw, [0], lw=0.5, label="")
    hline!(p3_hnw, [0], lw=0.5, label="")
    
    if is_gps
        
        println("GPS buoy")
        flush(stdout)
        
        # Locate GPS errors
        gps_errors = findall(isodd,parse.(Int,SubString.(string.(df.Column4), 2, 2), base = 16))

        if length(gps_errors) > 0
            vline!(p1_hnw, [wse_df.Date[1]], lw=0.5, ls=:dash, c=:red, label="GPS error")
            println(length(gps_errors)," GPS errors detected")
            flush(stdout)
        end

        for i in gps_errors
            vline!(p1_hnw, [wse_df.Date[i]], lw=0.5, ls=:dash, c=:red, label="")
        end
    
    end

    # get plotting limits
    x_lim1 = xlims(p1_hnw)[1]; y_lim1 = ylims(p1_hnw)[1]
    x_lim2 = xlims(p1_hnw)[2]; y_lim2 = ylims(p1_hnw)[2]

    # display plots to screen
    tm_tick = range(first(wse_df.Date),last(wse_df.Date),step=Minute(5))
    ticks = Dates.format.(tm_tick,"MM:SS")


    # display plots to screen
    plot_wse = Plots.plot(p1_hnw, p2_hnw, p3_hnw, layout = (3, 1), size = (1400, 600),
        xlim=(first(wse_df.Date),last(wse_df.Date)), xticks=(tm_tick,ticks), xtickfontsize=7,ytickfontsize=8,
        framestyle = :box,fg_legend=:transparent, legend=:bottomleft,
        leftmargin = 15Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1)            

    display(plot_wse)
    
    return()
    
    end    # plot_wses)


function do_fft(heave, N)
################################################
# calculate the Fourier coefficients vide (5.6.2)

    return([sum([heave[k]*exp(2*pi*-1im*k*l/N) for k in (1:N)]) for l in (1:N)])

    end    # do_fft()


function calc_psd(Hl, N)
################################################
# The power spectral density is obtained from the Fourier coefficients
    
    PSD = zeros(trunc(Int,N/2))

    for l = 1:trunc(Int,N/2)   
        if (l==1) || (l==trunc(Int,N/2)-1)
            PSD[l] = abs(Hl[l])^2
        else
            PSD[l] = abs(Hl[l])^2+abs(Hl[N-l-1])^2
        end
    end

    # Smooth coefficients vide (5.6.6)
    PSD_smooth = PSD
    [PSD_smooth[i] = PSD[i-1]/4 + PSD[i]/2 + PSD[i+1]/4 for i in (2:trunc(Int,N/2)-1)]

    return(PSD_smooth)

    end    # calc_psd()


function calc_tp5(f2,Sf)
##########################################
# Calculate Tp5 via Read method
    
    Sf_max = maximum(Sf)

    numerator = 0; denominator = 0

    Sf_sum = cumsum(Sf.*Sf_max).^5

    for i in eachindex(f2)
        w = Sf[i] / Sf_max
        numerator +=  f2[i] * w^5
        denominator += w^5
    end

    Fp5 = numerator / denominator
    
    return(Fp5)    # calc_tp5()

    end    # calc_tp5()


function calc_hm0(Sf,freq)
########################################## 
    
    ax1 = (last(freq) - first(freq)) / (length(freq)-1)

    # calc spectral moments m0, m1, m2, m3, and m4
    s00 = 0; m0 = 0

    for ii in 1:128

        s00 += freq[ii]^0 * Sf[ii];

    end

    m0 = 0.5*ax1*(first(freq)^0*first(Sf) + 2*s00 + last(freq)^0*last(Sf))

    return(4 * m0^0.5)

    end    # calc_hm0()


function calculate_frequency_domain_parameters(f2, spectra)
##########################################
# Calculate frequency-domain parameters    
# Calls: calc_tp5()
    
    ax1 = (last(f2) - first(f2)) / (length(f2)-1)

    # calc spectral moments m0, m1, m2, m3, and m4
    s00 = 0; s01 = 0; s02 = 0; s03 = 0; s04 = 0;
    m0 = 0; m1 = 0; m2 = 0; m3 = 0; m4 = 0

    for ii in 1:128

        s00 += f2[ii]^0 * spectra[ii]
        s01 += f2[ii]^1 * spectra[ii]
        s02 += f2[ii]^2 * spectra[ii]
        s03 += f2[ii]^3 * spectra[ii]
        s04 += f2[ii]^4 * spectra[ii]

    end

    m0 = 0.5*ax1*(first(f2)^0*first(spectra) + 2*s00 + last(f2)^0*last(spectra))
    m1 = 0.5*ax1*(first(f2)^1*first(spectra) + 2*s01 + last(f2)^1*last(spectra))
    m2 = 0.5*ax1*(first(f2)^2*first(spectra) + 2*s02 + last(f2)^2*last(spectra))
    m3 = 0.5*ax1*(first(f2)^3*first(spectra) + 2*s03 + last(f2)^3*last(spectra))
    m4 = 0.5*ax1*(first(f2)^4*first(spectra) + 2*s04 + last(f2)^4*last(spectra))

    ##println("m0 = ",m0," m1 = ",m1, " m2 = ",m2, " m3 = ",m2, " m4 = ",m4)

    # calc wave parameters Hm0, Hrms, T01, T02, Tc
    Hm0 = 4*sqrt(m0)     # Tucker & Pitt p.32 (2.2-6b)
    Hrms = sqrt(8*m0)    # Goda 2nd. Edition p.262 (9.15)
    T01 = m0/m1          # Tucker & Pitt p.41 Table 2.2 
    T02 = sqrt(m0/m2)    # Tucker & Pitt p.40 (2.3-2)
    Tc = sqrt(m2/m4)     # Tucker & Pitt p.41 Table 2.2 - also see Notes

    # identify spectral peak and frequency as peak
    Fp = f2[argmax(spectra)]
    Tp = 1/Fp
    fp5 = calc_tp5(f2, spectra)
    Tp5 = 1/fp5

    # calculate spectral width vide Tucker and Pitt p.85 (5.2-8)
    # Note: for JONSWAP, v = 0.39; for PM, v = 0.425
    v = (m0*m2 / m1^2 - 1)^0.5

    # calculate Skewness vide Tucker and Pitt p.109 (5.5-17)
    Skewness = (m0^2 * m3/m1^3 - 3*v^2 - 1) / v^3;
    
    return(Hm0, Hrms, T01, T02, Tc, Tp, fp5, Tp5, Skewness)
    
    end    # calculate_frequency_domain_parameters()


function calc_representative_spectra(frequency,Hm0,Tp,gamma)
##########################################    
    """
    function to calculate representative spectrum based on the Jonswap formula in Tucker and Pitt p.339 (10.3-9a)

    inputs:
        frequency - array of spectral frequencies
        Hm0 - floating point value
        Tp - floating point value
        gamma - floating point value - Peak ehhancement factor (Enter 1 for PM, or 3.3 for Jonswap)

        Typical calls:
        Spectra_PM = calc_representative_spectra(f2, Hm0, Tp, 1.0)
        Spectra_JONSWAP = calc_representative_spectra(f2, Hm0, Tp, 3.3)

    returns:
        Sf - array of representative spectra        
    """

    alpha = 1    # initial Philips constant (will decrease for each iteration required)
    g = 9.81
    fp = 1/Tp    # peak frequency

    hm0 = 99.    # set this to large value (so it will change on first iteration)

    Sf = [];

    while((Hm0 - hm0) <= 0.0005)

        Sf = vcat([alpha*g^2 * (2*pi)^-4 * ff^-5 * exp(-1.25 * (ff/fp)^-4) * gamma^exp(-(ff-fp)^2/(2*0.07^2 * fp^2)) for ff in frequency[findall(<=(fp), frequency)]],
                [alpha*g^2 * (2*pi)^-4 * ff^-5 * exp(-1.25 * (ff/fp)^-4) * gamma^exp(-(ff-fp)^2/(2*0.09^2 * fp^2)) for ff in frequency[findall(>(fp), frequency)]]);
        Sf[1] = 0;

###################################################################################################################################            
###  See discussion at https://stackoverflow.com/questions/44915116/how-to-decide-between-scipy-integrate-simps-or-numpy-trapz  ###
###################################################################################################################################

        ##        hm0 = 4*(np.trapz(Sf, frequency))^0.5    # calculate new Hm0 based on Sf values
        hm0 = calc_hm0(Sf,frequency);   # see https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.integrate.simps.html        
        alpha *= 0.95;    # reduce alpha by 5% so iterations approach a limit of 0.0005

    end

    return(Sf)
        
    end    # calc_representative_spectra()


function plot_spectra(wse_df)
################################################
    
    heave = wse_df.Heave
    Sample_frequency = 1.28
    
    # convert heave to matrix of individual 256-value spectra
    segments = arraysplit(heave, 512, 256)
    combined_segments = []
    
    for i in eachindex(segments)
        push!(combined_segments,power(periodogram(segments[i],nfft=512,fs=Sample_frequency,window=hanning)))
    end
    
    global freqs1 = freq(periodogram(segments[1],nfft=512,fs=Sample_frequency,window=hanning))
    global Pden = mean(combined_segments, dims = 1)

    # use Welch's method as a check
    global ps_w = welch_pgram(heave, 512, 256; onesided=true, nfft=512, fs=Sample_frequency, window=hanning);
    global f2 = freq(ps_w);
    global Pden2 = power(ps_w);

   Hm0, Hrms, T01, T02, Tc, Tp, fp5, Tp5, Skewness = calculate_frequency_domain_parameters(f2, Pden2)
    @printf("%s; Hm0 = %5.2fm; Hrms = %5.2fm; T01 = %5.2fs; T02 = %5.2fs; Tc = %5.2fs; Tp = %5.2fs; Tp5 = %5.2fs; Skewness = %5.4f",
        Dates.format(first(wse_df.Date), "yyyy-mm-dd HH:MM"),Hm0, Hrms, T01, T02, Tc, Tp, Tp5, Skewness)
    
    # Calculate representative spectra for P-M and JONSWAP
    Spectra_PM = calc_representative_spectra(f2, Hm0, Tp, 1.0);
    Spectra_JONSWAP = calc_representative_spectra(f2, Hm0, Tp, 3.3);

    # determing maximum y-axis value for spectral plots
    max_y = maximum([maximum(Pden[1]),maximum(Pden2),maximum(Spectra_JONSWAP)]) * 1.05
    
    if max_y < 0.1
        tick_val = 0.01
    elseif  max_y < 1
        tick_val = 0.1
    elseif max_y < 10
        tick_val = 1   
    else
        tick_val = 5
    end

    # Plot the representative spectra
    p_spectra = plot(f2,Spectra_JONSWAP, lw=2, c=:lightblue, label="JONSWAP spectrum (" * L"\gamma" * " = 3.3)")
    p_spectra = plot!(f2,Spectra_PM, lw=2, c=:lightgreen, label="Pierson-Moskowitz spectrum (" * L"\gamma" * " = 1.0)\n")
    
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    # Add frequency-domain parameters to plot
    x_lim = xlims(p_spectra)[1]; y_lim = ylims(p_spectra)[2]
    p_spectra = annotate!(x_lim*-25, y_lim*0.65, "Hm0 = " * string(round(Hm0, digits=2)) * "m",annotationfontsize=10) 
    p_spectra = annotate!(x_lim*-25, y_lim*0.60, "Hrms = " * string(round(Hrms, digits=2)) * "m") 
    p_spectra = annotate!(x_lim*-25, y_lim*0.55, "T01 = " * string(round(T01, digits=2)) * "s") 
    p_spectra = annotate!(x_lim*-25, y_lim*0.50, "T02 = " * string(round(T02, digits=2)) * "s") 
    p_spectra = annotate!(x_lim*-25, y_lim*0.45, "Tp = " * string(round(Tp, digits=2)) * "s") 
    p_spectra = annotate!(x_lim*-25, y_lim*0.40, "Tp5 = " * string(round(Tp5, digits=2)) * "s") 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
      
    # plot calculated spectra
    p_spectra = plot!(freqs1, Pden, label="Calc\n", 
        c=:yellow, lw=3, fillrange = 0, fillalpha = 0.05, fillcolor = :blue)
        
    # plot Welch's spectra
    p_spectra = plot!(f2, Pden2, label="Welch's method", 
        c=:red, lw=1, fillrange = 0, fillalpha = 0.05, fillcolor = :red)
    
    p_spectra = vline!([fp5; fp5], lw=1, ls =:dash, c=:red, label="Tp5")
    
    plot_spc = Plots.plot(p_spectra, layout = (1, 1), size = (1400, 600), framestyle = :box, 
        xlim=(0,0.64),  xticks = 0:0.05:1.28, xtickfontsize=7, ytickfontsize=8, xlabel="Frequency (Hertz)",
        ylim=(0,max_y), yticks=0:tick_val:max_y, ylabel="Spectral Density (sq.m/Hertz)",
        fg_legend=:transparent, title = " Spectral plot", titlefontsize=12,
        leftmargin = 15Plots.mm, bottommargin = 15Plots.mm, 
        grid=true, gridlinewidth=0.5, gridalpha=1, foreground_color_grid="lightgrey")            

    display(plot_spc)
    
    return()
    
    end    # plot_spectra()


function plot_spectrogram(wse_df)
    
    heave = wse_df.Heave;
    nw=256;
    spec = DSP.Periodograms.spectrogram(heave, nw, 250; fs=1.28,window=hanning);

    # display plots to screen
    tm_tick = range(first(wse_df.Date),last(wse_df.Date),step=Minute(5))
    ticks = Dates.format.(tm_tick,"MM:SS")

    spec1 = plot(first(wse_df.Date) + Microsecond.(ceil.((spec.time) * 1000000)), spec.freq, DSP.Periodograms.power(spec), lw=1, c=cgrad(:Spectral, rev=true), colorbar=false, 
        size=(1400, 600), framestyle = :box, title="Spectrogram", 
        xlim=(first(wse_df.Date),last(wse_df.Date)), xticks=(tm_tick,ticks), xtickfontsize=7, xlabel="Time (s)",
        ytickfontsize=8, ylabel="Frequency (Hz)",
        leftmargin = 15Plots.mm, bottommargin = 15Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1, show=true) 

    display(spec1)
    
    return()
    
    end    # plot_spectrogram()


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

    plot_wses(wse_df,is_gps)
    plot_spectrogram(wse_df)
    plot_spectra(wse_df)

end