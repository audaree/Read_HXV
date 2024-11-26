## Julia program to read a selected .HXV file and display 30-minute time series plots
## JW December 2022
#using ContinuousWavelets 
using CSV
using Dates, DataFrames, Distributions, DSP
using Gtk
using LaTeXStrings
using NativeFileDialog
using Printf
using Statistics #, StatsPlots
using Tk

function get_cyclic_word(word)
    ##########################################    
    # get cyclic data block - see 5.7.1.2 Spectrum fill or full wave spectrum
        
        a = ""
            for i in 1:4
                a = a * bitstring(parse(Int8, SubString.(word, i, i), base=16))[5:8]
            end
        
        return(a)
        
    end    # get_cyclic_word()
    
    
function process_system_file_word_0(system_file_word)
    ##########################################    
        
        Tp = parse(Int16, system_file_word[1:4]; base=2)
        M = parse(Int8, system_file_word[5]; base=2)
        T = parse(Int8, system_file_word[6]; base=2)
        F = parse(Int8, system_file_word[7]; base=2)
        C = parse(Int8, system_file_word[8]; base=2)
        Tn = parse(Int16, system_file_word[9:12]; base=2)
        
        println("Transmission No. ",Tn)
        println("System file word 0: Tp = ",Tp,"; M = ",M,"; T = ",T,"; F = ",F,"; C = ",C,"; Tn = ",Tn)
        
        return()
        
        end    # process_system_file_word_0()
             
        
function process_system_file_word_1(system_file_word)
##########################################    

    Hrms = parse(Int16, system_file_word[1:12]; base=2) / 400
    m0 = Hrms ^ 2
    Hm0 = 4 * sqrt(m0)
    @printf("System file word 1: Hrms = %5.2fm; m0 = %5.5f; Hm0 = %5.2fm\n",Hrms,m0,Hm0)
    push!(hm0,Hm0)

    return()

    end    # process_system_file_word_1()


function process_system_file_word_2(system_file_word)
##########################################    
    
    fz = parse(Int16, system_file_word[5:12]; base=2) / 400
    Tz = 1/fz
    @printf("System file word 2: fz = %5.4fHz; Tz = %5.2fs\n",fz,Tz)
    push!(tz,Tz)
    
    return()
    
    end    # process_system_file_word_2()


function process_system_file_word_3(system_file_word)
##########################################    
    
    global PSDmax = 5000 * exp(-parse(Int32, system_file_word[1:12]; base=2) / 200)
    @printf("System file word 3: PSDmax = %5.4fm²/Hz\n",PSDmax)
    push!(max_PSD,PSDmax)

    return()
    
    end    # process_system_file_word_3()


function process_system_file_word_4(system_file_word)
##########################################    
    
    Tr = parse(Int16, system_file_word[3:12]; base=2) / 20 - 5
    println("System file word 4: Tr = ",Tr,"ᵒC")

    return()
    
    end    # process_system_file_word_4()


function process_system_file_word_5(system_file_word)
##########################################    
    
    Tw = parse(Int16, system_file_word[3:12]; base=2) / 20 - 5
    println("System file word 5: Tw = ",Tw,"ᵒC")
    push!(tw,Tw)
    
    return()
    
    end    # process_system_file_word_5()


function process_system_file_word_6(system_file_word)
##########################################    
    
    B = parse(Int16, system_file_word[10:12]; base=2)
    tol = parse(Int16, system_file_word[1:8]; base=2)
    println("System file word 6: B = ",B,"; Tol = ",tol)

    return()
    
    end    # process_system_file_word_6()


function process_system_file_word_7(system_file_word)
##########################################    
    
    Av0 = parse(Int16, system_file_word[2:12]; base=2) / 800
    sign = parse(Int16, system_file_word[1]; base=2)
    
    if (sign==0)   
        println("System file word 7: Av0 = ",Av0,"m/s²")
    elseif (sign==1)
        println("System file word 7: Av0 = ",-Av0,"m/s²")
    end
    
    return()
    
    end    # process_system_file_word_7()


function process_system_file_word_8(system_file_word)
##########################################    
    
    Ax0 = parse(Int16, system_file_word[2:12]; base=2)
    sign = parse(Int16, system_file_word[1]; base=2)
    
    if (sign==0)   
        println("System file word 8: Ax0 = ",Ax0,"m/s²")
    elseif (sign==1)
        println("System file word 8: Ax0 = ",-Ax0,"m/s²")
    end
    
    return()
    
    end    # process_system_file_word_8()


function process_system_file_word_9(system_file_word)
##########################################    
    
    Ay0 = parse(Int16, system_file_word[2:12]; base=2)
    sign = parse(Int16, system_file_word[1]; base=2)
    
    if (sign==0)   
        println("System file word 9: Ay0 = ",Ay0,"m/s²")
    elseif (sign==1)
        println("System file word 9: Ay0 = ",-Ay0,"m/s²")
    end
    
    return()
    
    end    # process_system_file_word_9()


function process_system_file_word_10(system_file_word)
##########################################    
    
    LatMSB = system_file_word[2:12]
    SignLat = system_file_word[1]
    
##    println("System file word 10 done")
    
    return(LatMSB,SignLat)
    
    end    # process_system_file_word_10()


function process_system_file_word_11(system_file_word,LatMSB,SignLat)
##########################################    
    
    global LatLSB = system_file_word[1:12]
    
    lat = 90 * (parse(Int32, (LatMSB*LatLSB); base=2) / 2^23)
    if (SignLat=='0')   
        @printf("System file words 10 and 11: Latitude = %5.4fᵒN\n",lat)
        push!(latitude,lat)
    elseif (SignLat=='1')
        @printf("System file words 10 and 11: Latitude = %5.4fᵒS\n",lat)
        push!(latitude,-lat)        
    end
    
    return()
    
    end    # process_system_file_word_11()


function process_system_file_word_12(system_file_word)
##########################################    
    
    SignLon = system_file_word[1]
    LonMSB = system_file_word[2:12]
    
##    println("System file word 12 done")
    
    return(LonMSB,SignLon)
    
    end    # process_system_file_word_12()


function process_system_file_word_13(system_file_word,LonMSB,Sign)
##########################################    
    
    global LonLSB = system_file_word[1:12]
    
    lon = 180 * (parse(Int64, (LonMSB*LonLSB); base=2) / 2^23)
    
    if (Sign=='0')   
        @printf("System file words 12 and 13: Longitude = %5.4fᵒE\n",lon)
        push!(longitude,lon)
        
    elseif (Sign=='1')
        @printf("System file words 12 and 13: Longitude = %5.4fᵒW\n",lon)
        push!(longitude,-lon)
    end
    
    return()
    
    end    # process_system_file_word_13()


function process_system_file_word_14(system_file_word)
##########################################    
    
    O = 360 * (parse(Int16, system_file_word[5:12]; base=2) / 256)
    
    println("System file word 14: O = ",O,"ᵒ")
    
    return()
    
    end    # process_system_file_word_14()


function process_system_file_word_15(system_file_word)
##########################################    
    
    IncMSB = parse(Int16, system_file_word[5:12]; base=2)
    IncLSB = parse(Int16, system_file_word[1:4]; base=2)
    
    I = (90/128) * (IncMSB - 128 + IncLSB/16)
#    println("System file word 15: I = ",I,"ᵒ\n")
    @printf("System file word 15: I = %5.4fᵒ\n\n",I)
    
    return()
    
    end    # process_system_file_word_12()


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
    println("\nSelected ",infil,"\n")

    global df = DataFrame(CSV.File(infil,header=0, delim=",", types=String));

    system_df = DataFrame([name => Any[] for name in ["Word_No", "System_file_word"]])

    is_gps = false
    sync_word_location = findall(x -> x == "7FFF", df.Column2)

    for j in sync_word_location
        
        cyclic_word_2 = get_cyclic_word(df.Column2[j+1])
        Word_No = parse(Int, cyclic_word_2[1:4]; base=2)
        System_file_word = cyclic_word_2[5:16]
        push!(system_df,[Word_No,System_file_word])

    end

    global hm0 = []
    global tz = []
    global max_PSD = []
    global tw = []
    global latitude = []
    global longitude = []

    for i in eachindex(system_df.Word_No)
        if (system_df.Word_No[i]==0)
            process_system_file_word_0(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==1)
            process_system_file_word_1(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==2)
            process_system_file_word_2(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==3)
            process_system_file_word_3(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==4)
            process_system_file_word_4(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==5)
            process_system_file_word_5(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==6)
            process_system_file_word_6(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==7)
            process_system_file_word_7(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==8)
            process_system_file_word_8(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==9)
            process_system_file_word_9(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==10)
            global LatMSB,SignLat = process_system_file_word_10(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==11)
            try
                process_system_file_word_11(system_df.System_file_word[i],LatMSB,SignLat)
            catch
                println("Attempted to process System file word 11 - but no word 10 available")
            end
        elseif (system_df.Word_No[i]==12)
            global LonMSB,SignLon = process_system_file_word_12(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==13)
            try
                process_system_file_word_13(system_df.System_file_word[i],LonMSB,SignLon)
            catch
                println("Attempted to process System file word 13 - but no word 12 available")
            end
        elseif (system_df.Word_No[i]==14)
            process_system_file_word_14(system_df.System_file_word[i])
        elseif (system_df.Word_No[i]==15)
            process_system_file_word_15(system_df.System_file_word[i])
        else
            println("Error - Invalid System file word")       
        end
        
    end

    @printf("\nHm0 = %5.2fm; Tz = %5.2fs; Tw = %5.2fᵒC; Pdenₘₐₓ = %5.4fm²/Hz; Latitude = %5.4fᵒ; Longitude = %5.4fᵒ\n",
        mode(hm0),mode(tz),mode(tw),mode(max_PSD),mode(latitude),mode(longitude))
    println("-----------------------------------------------------------------------------------------------------------n\n")

end 