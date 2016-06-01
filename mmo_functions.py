"""Functions to model MMO engineering."""

import tellurium as te
import matplotlib.pyplot as plt

def set_sMMO_activation(sMMO_foldactivation, hours, *test_model):
    """Allows setting an sMMO fold-activation value and simulate methanol
    production, bioremediation and number of cells."""

    if test_model :
        model = test_model
        print "Test model loaded."
    else:
        model = '''
            J1: CH4 -> CH3OH; pMMO*v_pMMO + sMMO*v_sMMO
            J2: CH4 -> ; CH4*0.5/u
            J3: CH3OH -> ; CH3OH*0.5/u
            J4: VC -> bioremediation; n * (pMMO*v_pMMO_biorem + sMMO*v_sMMO_biorem)
            J5: n -> 2n; u # number of cells, need to implement stationary phase
        
            v_pMMO = Vmax_pMMO*CH4/(Km_pMMO + VC*Km_pMMO/Ki_pMMO + CH4)
            v_sMMO = Vmax_sMMO*CH4/(Km_sMMO + VC*Km_sMMO/Ki_sMMO + CH4)
        
            v_pMMO_biorem = Vmax_pMMO_VC*VC/(Ki_pMMO + CH4*Ki_pMMO/Km_pMMO + VC)
            v_sMMO_biorem = Vmax_sMMO_VC*VC/(Ki_sMMO + CH4*Ki_sMMO/Km_sMMO + VC)
        
            ###### Species init #######
            $CH4 = 0.00097 # 58uM/min; this is the methane uptake of 5GB1 cells
            CH3OH = 0 # arbitrarily set to zero
        
            ##### Param init ######
        
            # The following rates are for M. trichosporium  (Lee et al 2006)
            Vmax_pMMO = 0.082/60
            Km_pMMO = 8.3
            Vmax_pMMO_VC = 0.042/60
            Ki_pMMO = 26
        
            Vmax_sMMO = 0.726/60
            Km_sMMO = 92
            Vmax_sMMO_VC = 2.1/60
            Ki_sMMO = 160
        
            u = 6.5*3600 * k  # generation time * growth inhibition by compound
        
            pMMO = 1
            sMMO = 0.002 #500x repression in Cu medium with respect to no Cu medium
            $VC = 100
            k = pMMO * 0.69 + sMMO * 0.29
            n = 1 '''
        print "Using default model."

    models = []
    
    if test_model :     
        for i in sMMO_foldactivation:
            models.append(
                te.loada(
                    str.replace(test_model, "sMMO = 0.002", "sMMO = 0.002 * " + str(i))
                        )
                )
    else:
        for i in sMMO_foldactivation:
            models.append(
                te.loada(
                    str.replace(model, "sMMO = 0.002", "sMMO = 0.002 * " + str(i))
                        )
                )
        
    results = []
    for i in models:
        results.append(i.simulate(0, 3600*hours, 1000))

    time = []
    for i in results[0][:, 0]:
        time.append(i/3600)

    ylabels = ['CH3OH [uM]', 'Cumulative culture degradation [uM]','# cells']
    linestyles = ['r--', 'b--', 'g--', 'c--', 'm--', 'y--', 'k--']
    for i in range(3): # number of species to plot, set labels accordingly above
        sim = 0
        for result in results:
            plt.plot(time, result[:, i+1], linestyles[sim], label=str(sMMO_foldactivation[sim]))
            sim +=1
        plt.legend(loc='upper left', title = "sMMO fold-activation")
        plt.ylabel(ylabels[i])
        plt.xlabel('time [h]')
        plt.show()

def set_pMMO(pMMO, hours) :

    model = '''
    J1: CH4 -> CH3OH; pMMO*v_pMMO + sMMO*v_sMMO
    J2: CH4 -> ; CH4*0.5/u
    J3: CH3OH -> ; CH3OH*0.5/u
    J4: VC -> bioremediation; n * (pMMO*v_pMMO_biorem + sMMO*v_sMMO_biorem)
    J5: n -> 2n; u #number of cells

    # Assuming Ki is the Km for the inhibitor.
    # This does not look at bioremediation rates, just biomass production.

    v_pMMO = Vmax_pMMO*CH4/(Km_pMMO + VC*Km_pMMO/Ki_pMMO + CH4)
    v_sMMO = Vmax_sMMO*CH4/(Km_sMMO + VC*Km_sMMO/Ki_sMMO + CH4)

    v_pMMO_biorem = Vmax_pMMO_VC*VC/(Ki_pMMO + CH4*Ki_pMMO/Km_pMMO + VC)
    v_sMMO_biorem = Vmax_sMMO_VC*VC/(Ki_sMMO + CH4*Ki_sMMO/Km_sMMO + VC)

    ###### Species init #######
    $CH4 = 0.00097 #58uM/min
    CH3OH = 0

    ##### Param init ######
    Vmax_pMMO = 0.082/60
    Km_pMMO = 8.3
    Vmax_pMMO_VC = 0.042/60
    Ki_pMMO = 26

    Vmax_sMMO = 0.726/60
    Km_sMMO = 92
    Vmax_sMMO_VC = 2.1/60
    Ki_sMMO = 160

    u = 6.5*3600 * k  # generation time * growth inhibition by compound

    pMMO = 1
    sMMO = 1 - pMMO
    $VC = 10
    k = pMMO * 0.69 + sMMO * 0.29
    n = 1

'''

    models = []
    for i in pMMO:
        models.append(
            te.loada(
                str.replace(model, "pMMO = 1", "pMMO =" + str(i) )
                    )
            )


    results = []
    for i in models:
        results.append(i.simulate(0, 3600*hours, 1000))

    time =[]
    for i in results[0][:,0]: time.append(i/3600)


    ylabels = ['CH3OH [uM]', 'Cumulative culture degradation [uM]', '# of cells']
    linestyles = ['r--', 'b--', 'g--', 'c--', 'm--', 'y--', 'k--']
    for i in range(3): #number of species to plot, set labels accordingly above
        sim = 0
        for result in results:
            plt.plot(time, result[:,i+1], linestyles[sim], label=str(pMMO[sim]))
            sim +=1
        plt.legend(loc='upper left', title = "pMMO %")
        plt.ylabel(ylabels[i])
        plt.xlabel('time [h]')
        plt.show()


def all_pMMO_sMMO(hours) :
    set_pMMO([0,1], hours)
