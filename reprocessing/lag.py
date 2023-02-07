import PYCCF_v1 as myccf
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os


"""
judge time lags is calculated??
"""
def calculate_lags(filepath, w_band, plot=False):
    lcs = fits.open(filepath)

    """
    plot light curves.
    """
    if plot == True:
        plt.figure(figsize=(15, 10))

        colors= ['black', 'red', 'orange', 'green', 'cyan', 'purple']

        for i in range(len(w_band)):
            mjd1 = lcs[1].data['time'].flatten()
            flux1 = lcs[1].data[w_band[i]].flatten()

            plt.plot(mjd1, flux1, lw=2, label=w_band[i], color=colors[i])

        plt.xlabel('MJD (days)', fontsize=20)
        plt.ylabel('Flux', fontsize=20)
        plt.legend(ncol=3)
        plt.grid()
        plt.tick_params(top='on', right='on', which='major',labelsize=20, width=2, length=2)
        plt.tick_params(top='on', right='on', which='minor',labelsize=20, width=2, length=5)
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        # plt.savefig("../dat/six_band_lc.png", dpi=300)
        # plt.show()
        # plt.close()

    """
    calculate the time lags.
    """
    
    sim_lags_swuw2 = []
    ccf = []
    taulist = []
    cols = []
    colccf = []
    coltau = []
    for i_band in range(len(w_band)):
        pack = myccf.peakcent(lcs[1].data['time'].flatten(), lcs[1].data[w_band[1]].flatten(),
                                lcs[1].data['time'].flatten(), lcs[1].data[w_band[i_band]].flatten(),
                                -50, 50, 0.01, thres=0.8, siglevel=0.95,
                                imode=0, sigmode=0.2)
        sim_lags_swuw2.append(round(pack[0],3))
        ccf.append(pack[4][0])
        taulist.append(pack[4][1])

    for i in range(len(w_band)):
        print(i)
        cols.append(
            fits.Column(name=w_band[i]+'_lags_swuw2', format='D', array=np.array(sim_lags_swuw2[0+i:1+i]))
            )
        colccf.append(
            fits.Column(name=w_band[i]+'_ccf', format='D', array=np.array(ccf[i]))
            )
        coltau.append(
            fits.Column(name=w_band[i]+'_tau', format='D', array=np.array(taulist[i]))
            )

    lcs.append(fits.BinTableHDU.from_columns(fits.ColDefs(cols),name="time lag"))
    lcs.append(fits.BinTableHDU.from_columns(fits.ColDefs(colccf), name="ccf-r"))
    lcs.append(fits.BinTableHDU.from_columns(fits.ColDefs(coltau),name="ccf-tau"))
    lcs.writeto(filepath, overwrite=True)
    lcs.close()
    print(" ")
    # print("Now, the time lags result have saved to: "+folder+"/simulationInfo.fits")

    return True



if __name__ == '__main__':
    import sys
    import configparser
    import ast

    # config = configparser.ConfigParser()
    # try:
    #     config.read('DiskSimulation.ini')
    # except IOError:
    #     sys.exit("The configuration file not found.")

    folder = str(sys.argv[1])
    # w_band = ast.literal_eval(config['miscellaneous']['w_band'])
    w_band = []
    # band = sys.argv[2:]
    # band = ["u", "g", "r", "i", "z"]
    band = ["uw2", "um2", "uw1", "uuu", "ubb", "uvv"]

    for i in range(len(band)):
        w_band.append(band[i])

    calculate_lags(filepath=folder, w_band=w_band, plot=False)
    