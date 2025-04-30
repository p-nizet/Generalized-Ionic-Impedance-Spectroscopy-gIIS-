def Wo(w_values, R, Ri, C, gIIS_type):
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    pi = math.pi
    w = w_values

    wRC = w*Ri*C
    x   = np.sqrt(1j*wRC)
    sinh= (np.exp(x)-np.exp(-1*x)) / 2
    cosh= (np.exp(x)+np.exp(-1*x)) / 2
    tanh= (np.exp(x)-np.exp(-x)) / (np.exp(x)+np.exp(-x))
    Warb= Ri / (tanh * x)
    
    if gIIS_type == 'total':
        print(f"✅ Executing Wo total gIIS with R_s = {R}, R_i = {Ri}, C_chem = {C}")

        QV = tanh / x

        H= Warb/(R+Warb)
        OIS= H*QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)


        # Crear la figura y el primer eje
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')

        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(0, 1)
        plt.ylim(0, 0.5)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'local':
        print(f"✅ Executing Wo local gIIS with R_s = {R}, R_i = {Ri}, C_chem = {C}")
        Position = float(input("Write the position (in x/L units) where you measure: "))

        coshkLx = np.cosh(x*(1-Position))
        QV = coshkLx / cosh

        H= Warb/(R+Warb)
        OIS= H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi + Ang[i]


        # Crearte figures
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 180)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.1, 1)
        plt.ylim(0, 0.6)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'partial':
        print(f"✅ Executing Wo partial gIIS with R_s = {R}, R_i = {Ri}, C_chem = {C}")
        Range1 = float(input("Write the initial position (in x/L units): "))
        Range2 = float(input("Write the final position (in x/L units): "))

        # Create a figure to ilustrate the area under study
        fig, ax = plt.subplots(figsize=(8, 1))
        ax.fill_betweenx([0.4, 0.6], Range1, Range2, color='skyblue', alpha=0.6)
        ax.set_xlim(0, 1)
        ax.set_ylim(0.4, 0.6)
        ax.set_xticks([])  # Hide x-axis ticks
        ax.set_yticks([])  # Hide y-axis ticks
        ax.set_title('Observed area of the material', fontsize=14)
        ax.set_xlabel('x/L', fontsize=12)
        plt.show()

        sinhA = np.sinh((1-Range2)*x)
        sinhB = np.sinh((1-Range1)*x)
        QV = (1/(Range2 - Range1)) / x * (sinhB-sinhA) / (cosh)

        H = Warb / (R + Warb)
        OIS = H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi/2 - Ang[i]



        # Create figures
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 180)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'eis':
        
        EIS= R + Warb

        # Extraer la parte real e imaginaria
        result_real = np.real(EIS)
        result_imag = np.imag(EIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for val in Ang:
            if val < 0:
                val = pi/2 + val
            else:
                val=val


        # Create Bode figure
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode EIS plot')
        plt.show()

        # Create Nyquist figure
        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist EIS plot')
        plt.xlim(0, R+Ri)
        plt.ylim(-1*(R+Ri)/20, (R+Ri)/2)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()
        








def Wo_el(w_values, R, Ri, Re, C, gIIS_type):
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    pi = math.pi
    w = w_values

    wRC = w*(Ri+Re)*C
    x   = np.sqrt(1j*wRC)
    sinh= (np.exp(x)-np.exp(-x)) / 2
    cosh= (np.exp(x)+np.exp(-x)) / 2
    tanh= sinh/cosh
    Warb= ((Ri**2+Re**2)*cosh+2*Re*Ri+Re*Ri*x*sinh)/(x*(Re+Ri)*sinh)

    
    if gIIS_type == 'total':
        print(f"✅ Executing Wo total gIIS with R_s = {R}, R_i = {Ri}, R_e = {Re}, C_chem = {C}")

        QV = (Ri+Re)/x*(Ri+Re)*(sinh)/((Ri**2+Re**2)*cosh+2*Re*Ri+Re*Ri*x*sinh)

        H= Warb/(R+Warb)
        OIS= H*QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for val in Ang:
            if val < 0:
                val = pi/2 + val
            else:
                val=val


        # Create Bode figure
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        # Create Nyquist figure
        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'local':
        print(f"✅ Executing Wo local gIIS with R_s = {R}, R_i = {Ri}, R_e = {Re}, C_chem = {C}")
        Position = float(input("Write the position (in x/L units) where you measure: "))

        coshkLx = np.cosh(x*(1-Position))
        coshkx  = np.cosh(x*Position)
        QV = (Ri+Re)*(Ri*coshkLx+Re*coshkx)/((Ri**2+Re**2)*cosh+2*Re*Ri+Re*Ri*x*sinh)

        H= Warb/(R+Warb)
        OIS= H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi + Ang[i]


        # Create figures
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 180)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.1, 1)
        plt.ylim(0, 0.6)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'partial':
        print(f"✅ Executing Wo partial gIIS with R_s = {R}, R_i = {Ri}, R_e = {Re}, C_chem = {C}")
        Range1 = float(input("Write the initial position (in x/L units): "))
        Range2 = float(input("Write the final position (in x/L units): "))

        # Create a figure to ilustrate the area under study
        fig, ax = plt.subplots(figsize=(8, 1))
        ax.fill_betweenx([0.4, 0.6], Range1, Range2, color='skyblue', alpha=0.6)
        ax.set_xlim(0, 1)
        ax.set_ylim(0.4, 0.6)
        ax.set_xticks([])  # Hide x-axis ticks
        ax.set_yticks([])  # Hide y-axis ticks
        ax.set_title('Observed area of the material', fontsize=14)
        ax.set_xlabel('x/L', fontsize=12)
        plt.show()

        sinha   = np.sinh(x*(1-Range1))
        sinhb   = np.sinh(x*(1-Range2))
        sinhA   = np.sinh(x*Range1)
        sinhB   = np.sinh(x*Range2)

        QV = (1/(Range2 - Range1)) / x * (Ri+Re)*(Ri*(sinha-sinhb)+Re*(sinhB-sinhA))/((Ri**2+Re**2)*cosh+2*Re*Ri+Re*Ri*x*sinh)


        H = Warb / (R + Warb)
        OIS = H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi/2 - Ang[i]



        # Crear la figura y el primer eje
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 180)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'eis':
        
        EIS= R + Warb

        # Extraer la parte real e imaginaria
        result_real = np.real(EIS)
        result_imag = np.imag(EIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for val in Ang:
            if val < 0:
                val = pi/2 + val
            else:
                val=val


        # Create Bode figure
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode EIS plot')
        plt.show()

        # Create Nyquist figure
        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist EIS plot')
        plt.xlim(0, R+Ri+Re)
        plt.ylim(-1*(R+Ri+Re)/20, (R+Ri+Re)/2)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()




def Wo_pore(w_values, R, Ri, Re, C, RiW, Rinc, Cdl, gIIS_type):
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    pi = math.pi
    w = w_values

    #Particle warburg
    wRC = w*RiW*C
    x   = np.sqrt(1j*wRC)
    sinh= (np.exp(x)-np.exp(-1*x)) / 2
    cosh= (np.exp(x)+np.exp(-1*x)) / 2
    tanh= sinh/cosh
    W   = RiW / (tanh * x)

    #Total transmission line
    ZC_sust= (1j*w*Cdl + 1/(Rinc + W))**(-1)
    wRC2 = (Ri+Re) / ZC_sust
    x2   = np.sqrt(wRC2)
    sinh2= np.sinh(x2)
    cosh2= np.cosh(x2)
    tanh2 = sinh2 / cosh2
    Warb= ((Ri**2+Re**2)*cosh2+2*Re*Ri+Re*Ri*x2*sinh2)/(x2*(Re+Ri)*sinh2)

    
    if gIIS_type == 'total':
        print(f"✅ Executing Wo_pore total gIIS with R_s = {R}, R_pore = {Ri}, R_e = {Re}, C_chem = {C}, R_ipore = {RiW}, R_inc = {Rinc}, C_dl = {Cdl}")

        Q = (Ri+Re)/x2*(Ri+Re)*(sinh2)/((Ri**2+Re**2)*cosh2+2*Re*Ri+Re*Ri*x2*sinh2)

        QV= Q * W/(Rinc+W) * tanh / x

        H= Warb/(R+Warb)
        OIS= H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for val in Ang:
            if val < 0:
                val = pi/2 + val
            else:
                val=val


        # Create figures
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()



    elif gIIS_type == 'local':
        print(f"✅ Executing Wo_pore local gIIS with R_s = {R}, R_pore = {Ri}, R_e = {Re}, C_chem = {C}, R_ipore = {RiW}, R_inc = {Rinc}, C_dl = {Cdl}")
        Position = float(input("Write the position (in x/L units) where you measure: "))

        coshkLx = np.cosh(x2*(1-Position))
        coshkx  = np.cosh(x2*Position)
        Q = (Ri+Re) * (Ri*coshkLx + Re*coshkx) /((Re**2+Ri**2)*cosh2+2*Ri*Re+Ri*Re*x2*sinh2)

        QV= Q * W/(Rinc+W) * tanh / x

        H= Warb/(R+Warb)
        OIS= H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi + Ang[i]


        # Create figures
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 180)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()



    elif gIIS_type == 'partial':
        print(f"✅ Executing Wo_pore partial gIIS with R_s = {R}, R_pore = {Ri}, R_e = {Re}, C_chem = {C}, R_ipore = {RiW}, R_inc = {Rinc}, C_dl = {Cdl}")
        Range1 = float(input("Write the initial position (in x/L units): "))
        Range2 = float(input("Write the final position (in x/L units): "))

        # Create a figure to ilustrate the area under study
        fig, ax = plt.subplots(figsize=(8, 1))
        ax.fill_betweenx([0.4, 0.6], Range1, Range2, color='skyblue', alpha=0.6)
        ax.set_xlim(0, 1)
        ax.set_ylim(0.4, 0.6)
        ax.set_xticks([])  # Hide x-axis ticks
        ax.set_yticks([])  # Hide y-axis ticks
        ax.set_title('Observed area of the material', fontsize=14)
        ax.set_xlabel('x/L', fontsize=12)
        plt.show()

        sinha   = np.sinh(x2*(1-Range1))
        sinhb   = np.sinh(x2*(1-Range2))
        sinhA  = np.sinh(x2*Range1)
        sinhB  = np.sinh(x2*Range2)

        QV = (1/(Range2 - Range1)) / x2 * (Ri+Re)*(Ri*(sinha-sinhb)+Re*(sinhB-sinhA))/((Ri**2+Re**2)*cosh2+2*Re*Ri+Re*Ri*x2*sinh2)


        H = Warb / (R + Warb)
        OIS = H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi/2 - Ang[i]



        # Create figures
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 180)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'eis':
        
        EIS= R + Warb

        # Extraer la parte real e imaginaria
        result_real = np.real(EIS)
        result_imag = np.imag(EIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for val in Ang:
            if val < 0:
                val = pi/2 + val
            else:
                val=val


        # Create Bode figure
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode EIS plot')
        plt.show()

        # Create Nyquist figure
        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist EIS plot')
        plt.xlim(0, R+Rinc+Ri+Re+RiW/2)
        plt.ylim(-1*(R+Rinc+Ri+Re+RiW/2)/20, (R+Rinc+Ri+Re+RiW/2)/2)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()








def Ws(w_values, R, Ri, Re, Rinc, C, gIIS_type):
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    pi = math.pi
    w = w_values

    wRC = w*(Ri+Re)*C
    x   = np.sqrt(1j*wRC)
    sinh= (np.exp(x)-np.exp(-1*x)) / 2
    cosh= (np.exp(x)+np.exp(-1*x)) / 2
    tanh= sinh/cosh
    Warb= ((Ri**2+Re**2+1/Rinc*(Re*Ri**2+Ri*Re**2))*cosh+2*Re*Ri+(Re*Ri*x+Ri/Rinc*(Ri**2+Re*Ri)/x)*sinh) / (x*(Re+Ri)*sinh + (cosh*(Re+Ri)**2) / Rinc)
    
    if gIIS_type == 'total':
        print(f"✅ Executing Ws total gIIS with R_s = {R}, R_i = {Ri}, R_e = {Re}, R_inc = {Rinc}, C_chem = {C}")

        QV = (Re+Ri)**2 / x * (sinh + Ri/(Rinc*x)*(cosh - 1)) / ((Ri**2+Re**2+1/Rinc*(Re*Ri**2+Ri*Re**2))*cosh+2*Re*Ri+(Re*Ri*x+Ri/Rinc*(Ri**2+Re*Ri)/x)*sinh)

        H= Warb / (R + Warb)
        OIS= H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)


        # Crear la figura y el primer eje
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        #plt.xlim(-0.1, 1)
        #plt.ylim(-0.05, 0.5)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()
        
        


    elif gIIS_type == 'local':
        print(f"✅ Executing Ws local gIIS with R_s = {R}, R_i = {Ri}, R_e = {Re}, R_inc = {Rinc}, C_chem = {C}")
        Position = float(input("Write the position (in x/L units) where you measure: "))

        coshkLx = np.cosh(x*(1-Position))
        coshkx  = np.cosh(x*Position)
        sinhkLx = np.sinh(x*(1-Position))
        QV = (Ri+Re)*(Ri*coshkLx + Re*coshkx + Ri/(x*Rinc)*(Ri+Re)*sinhkLx)/((Ri**2+Re**2+1/Rinc*(Re*Ri**2+Ri*Re**2))*cosh+2*Re*Ri+(Re*Ri*x+Ri/Rinc*(Ri**2+Re*Ri)/x)*sinh)

        H= Warb/(R+Warb)
        OIS= H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi + Ang[i]


        # Crear la figura y el primer eje
        fig, ax1 = plt.subplots()

        # Graficar la primera serie de datos
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        # Segundo eje comparte el mismo eje x y grafico la segunda serie de datos
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 180)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.1, 1)
        plt.ylim(-0.05, 0.6)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()




    elif gIIS_type == 'partial':
        print(f"✅ Executing Ws partial gIIS with R_s = {R}, R_i = {Ri}, R_e = {Re}, R_inc = {Rinc}, C_chem = {C}")
        Range1 = float(input("Write the initial position (in x/L units): "))
        Range2 = float(input("Write the final position (in x/L units): "))

        # Create a figure to ilustrate the area under study
        fig, ax = plt.subplots(figsize=(8, 1))
        ax.fill_betweenx([0.4, 0.6], Range1, Range2, color='skyblue', alpha=0.6)
        ax.set_xlim(0, 1)
        ax.set_ylim(0.4, 0.6)
        ax.set_xticks([])  # Hide x-axis ticks
        ax.set_yticks([])  # Hide y-axis ticks
        ax.set_title('Observed area of the material', fontsize=14)
        ax.set_xlabel('x/L', fontsize=12)
        plt.show()

        sinha   = np.sinh(x*(1-Range1))
        sinhb   = np.sinh(x*(1-Range2))
        sinhA   = np.sinh(x*Range1)
        sinhB   = np.sinh(x*Range2)
        cosha   = np.cosh(x*(1-Range1))
        coshb   = np.cosh(x*(1-Range2))

        QV = 1/(Range2-Range1) / x  * (Ri+Re)*(Ri*(sinha-sinhb) + Re*(sinhB-sinhA) + Ri/(x*Rinc)*(Ri+Re)*(cosha-coshb))/((Ri**2+Re**2+1/Rinc*(Re*Ri**2+Ri*Re**2))*cosh+2*Re*Ri+(Re*Ri*x+Ri/Rinc*(Ri**2+Re*Ri)/x)*sinh)

        H= Warb/(R+Warb)
        OIS= H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi + Ang[i]


        # Create figure
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 180)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.1, 1)
        plt.ylim(-0.05, 0.6)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'eis':
        
        EIS= R + Warb

        # Extraer la parte real e imaginaria
        result_real = np.real(EIS)
        result_imag = np.imag(EIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for val in Ang:
            if val < 0:
                val = pi/2 + val
            else:
                val=val


        # Create Bode figure
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode EIS plot')
        plt.show()

        # Create Nyquist figure
        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist EIS plot')
        plt.xlim(0, R+Rinc+Ri+Re)
        plt.ylim(-1*(R+Rinc+Ri+Re)/20, (R+Rinc+Ri+Re)/2)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()
        






def SOFC(w_values, R, Ri, Rinc, C, gIIS_type):
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    pi = math.pi
    w = w_values

    iwRC = Ri*(1j*w*C + 1/Rinc)
    x   = np.sqrt(iwRC)
    sinh= (np.exp(x)-np.exp(-1*x)) / 2
    cosh= (np.exp(x)+np.exp(-1*x)) / 2
    tanh= sinh/cosh
    Warb= Ri / (tanh * x)

    
    if gIIS_type == 'total':
        print(f"✅ Executing SOFC electrode total gIIS with R_s = {R}, R_i = {Ri}, R_inc = {Rinc}, C_chem = {C}")

        QV = tanh / x

        H= Warb/(R+Warb)
        OIS= H*QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)

        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)


        # Crear la figura y el primer eje
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()
    


    elif gIIS_type == 'local':
        print(f"✅ Executing SOFC electrode local gIIS with R_s = {R}, R_i = {Ri}, R_inc = {Rinc}, C_chem = {C}")
        Position = float(input("Write the position (in x/L units) where you measure: "))

        coshkLx = np.cosh(x*(1-Position))
        QV = coshkLx / cosh

        H= Warb/(R+Warb)
        OIS= H * QV

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi + Ang[i]


        # Crear la figura y el primer eje
        fig, ax1 = plt.subplots()

        # Graficar la primera serie de datos
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        # Segundo eje comparte el mismo eje x y grafico la segunda serie de datos
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 180)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.1, 1)
        plt.ylim(0, 0.6)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()



    elif gIIS_type == 'partial':
        print(f"✅ Executing SOFC electrode partial gIIS with R_s = {R}, R_i = {Ri}, R_inc = {Rinc}, C_chem = {C}")
        Range1 = float(input("Write the initial position (in x/L units): "))
        Range2 = float(input("Write the final position (in x/L units): "))

        # Create a figure to ilustrate the area under study
        fig, ax = plt.subplots(figsize=(8, 1))
        ax.fill_betweenx([0.4, 0.6], Range1, Range2, color='skyblue', alpha=0.6)
        ax.set_xlim(0, 1)
        ax.set_ylim(0.4, 0.6)
        ax.set_xticks([])  # Hide x-axis ticks
        ax.set_yticks([])  # Hide y-axis ticks
        ax.set_title('Observed area of the material', fontsize=14)
        ax.set_xlabel('x/L', fontsize=12)
        plt.show()


        sinh34 = np.sinh((1-Range1)*x)
        sinh12 = np.sinh((1-Range2)*x)
        QV = (1/(Range2-Range1)) / x * (sinh34 - sinh12) / cosh

        H = Warb / (R + Warb)
        OIS = H * QV



        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)

        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi/2 - Ang[i]


        # Crear la figura y el primer eje
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'eis':
        
        EIS= R + Warb

        # Extraer la parte real e imaginaria
        result_real = np.real(EIS)
        result_imag = np.imag(EIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for val in Ang:
            if val < 0:
                val = pi/2 + val
            else:
                val=val


        # Create Bode figure
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode EIS plot')
        plt.show()

        # Create Nyquist figure
        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist EIS plot')
        plt.xlim(0, R+Rinc+Ri)
        plt.ylim(-1*(R+Rinc+Ri)/20, (R+Rinc+Ri)/2)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()









def Two_C(w_values, R, R1, C1, R2, C2, gIIS_type):
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    pi = math.pi
    w = w_values


    if gIIS_type == 'total':
        Vision = input("👉 Write what species you want to see (1, 2, or both): ").strip().lower()
        print(f"IMPORTANT⚠️: Executing Two species gIIS without diffusion with normalization to 1")
        
        if Vision == '1':
            f1 = 1
            f2 = 0
        elif Vision == '2':
            f1 = 0
            f2 = 1
        elif Vision == 'both':        
            f1 = C1/(C1+C2)
            f2 = C2/(C1+C2)

        Z1= 1/(1j*w*C1)
        Z2= 1/(1j*w*C2)

        # Transfer functions
        H1= Z1 / (R1+Z1+R*(1+(Z1+R1)/(Z2+R2)))
        H2= Z2 / (R2+Z2+R*(1+(Z2+R2)/(Z1+R1)))

        #Total gIIS
        OIS  = H1*f1 + H2*f2


        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)

        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi/2 - Ang[i]


        # Crear la figura y el primer eje
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode Total gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist Total gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'eis':
        
        EIS= R + (1/(R1+(1j*w*C1)**(-1))+1/(R2+(1j*w*C2)**(-1)))**(-1)

        # Extraer la parte real e imaginaria
        result_real = np.real(EIS)
        result_imag = np.imag(EIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for val in Ang:
            if val < 0:
                val = pi/2 + val
            else:
                val=val


        # Create Bode figure
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode EIS plot')
        plt.show()

        # Create Nyquist figure
        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist EIS plot')
        plt.xlim(0, R+(R1+R2+Ri1+Ri2)/2)
        plt.ylim(-1*(R+(R1+R2+Ri1+Ri2)/2)/20, (R+(R1+R2+Ri1+Ri2)/2)/2)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()













def Two_W(w_values, R, R1, Ri1, C1, R2, Ri2, C2, gIIS_type):
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    pi = math.pi
    w = w_values

    #Wo1
    wRC1 = w*Ri1*C1
    x1   = np.sqrt(1j*wRC1)
    tanh1= (np.exp(x1)-np.exp(-x1)) / (np.exp(x1)+np.exp(-x1))
    Warb1= Ri1 / (tanh1 *x1)
    QV1 = tanh1 / x1

    #Wo2
    wRC2 = w*Ri2*C2
    x2   = np.sqrt(1j*wRC2)
    tanh2= (np.exp(x2)-np.exp(-x2)) / (np.exp(x2)+np.exp(-x2))
    Warb2= Ri2 / (tanh2 *x2)
    QV2 = tanh2 / x2


    if gIIS_type == 'total':
        Vision = input("👉 Write what species you want to see (1, 2, or both): ").strip().lower()
        print(f"IMPORTANT⚠️: Executing Two species gIIS without diffusion with normalization to 1")
        if Vision == '1':
            f1 = 1
            f2 = 0
        elif Vision == '2':
            f1 = 0
            f2 = 1
        elif Vision == 'both':        
            f1 = C1/(C1+C2)
            f2 = C2/(C1+C2)

        # Transfer functions
        H1= Warb1 / (R1+Warb1+R*(1+(Warb1+R1)/(Warb2+R2)))
        H2= Warb2 / (R2+Warb2+R*(1+(Warb2+R2)/(Warb1+R1)))

        #Total gIIS
        OIS  = H1*f1*QV1 + H2*f2*QV2

        # Extraer la parte real e imaginaria
        result_real = np.real(OIS)
        result_imag = np.imag(OIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)

        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for i in range(len(Ang)):
            if Ang[i] < 0:
                Ang[i] = pi/2 - Ang[i]


        # Crear la figura y el primer eje
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_ylim(0, 1.05)
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode Total gIIS plot')
        plt.show()

        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist Total gIIS plot')
        plt.xlim(-0.15, 1)
        plt.ylim(-0.05, 0.55)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()


    elif gIIS_type == 'eis':
        
        EIS= R + (1/(R1+Warb1)+1/(R2+Warb2))**(-1)

        # Extraer la parte real e imaginaria
        result_real = np.real(EIS)
        result_imag = np.imag(EIS)

        # Imprimir los resultados
        ParteReal=[]
        ParteImag=[]
        for w, real, imag in zip(w_values, result_real, result_imag):
            ParteReal.append(real)
            ParteImag.append(-1*imag)


        ParteReal = np.array(ParteReal)
        ParteImag = np.array(ParteImag)
        Arg=(ParteReal**2+ParteImag**2)**0.5
        Ang=np.arctan(ParteImag/ParteReal)
        for val in Ang:
            if val < 0:
                val = pi/2 + val
            else:
                val=val


        # Create Bode figure
        fig, ax1 = plt.subplots()
        ax1.plot(w_values/(2*np.pi), Arg, '--r', label='Amp')
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax2 = ax1.twinx()
        ax2.plot(w_values/(2*np.pi), Ang*180/pi, '--g', label='Diphase')
        ax2.set_ylim(0, 95)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Amplitude')
        ax2.set_ylabel('-Diphase')
        fig.legend(loc='upper right')
        plt.title('Bode EIS plot')
        plt.show()

        # Create Nyquist figure
        plt.plot(ParteReal, ParteImag)
        plt.title('Nyquist EIS plot')
        plt.xlim(0, R+(R1+R2+Ri1+Ri2)/2)
        plt.ylim(-1*(R+(R1+R2+Ri1+Ri2)/2)/20, (R+(R1+R2+Ri1+Ri2)/2)/2)
        plt.axhline(y=0, color='grey', linewidth=0.5)
        plt.xlabel('Real')
        plt.ylabel('-Imaginary')
        plt.show()