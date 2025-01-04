#######################################################
#                 SafeGuard  Analysis                 #
#######################################################
# Labib Muhammad Majdi
# Zaki Adzani Sutrisno
# Thirdy Aji W. Pamungkas 
#######################################################
# Referensi utama
# Shugart, N., & King, J. (2017). A new modeling technique to analyze safeguards measurements in large systems. Nuclear Technology, 199(2), 129-150.

#######################################################
# import module
import numpy as np
import random
import matplotlib.pyplot as plt

#######################################################
# data SQ
SQ_LEU      = 75.
SQ_HEU      = 25.
SQ_233U     = 8.
SQ_238PU    = 80.

#######################################################
# class
class IKMP:
    # definition
    def __init__(self, mass: float, sigma_sd: float, sigma_so: float, sigma_r: float, t_cal: int = 0, drift: float = 0., offset: float = 0., rand: float = 0., type: str = 'flow'):
        # uncertainty
        self.sigma_sd   = sigma_sd
        self.sigma_so   = sigma_so
        self.sigma_r    = sigma_r
        # mass
        self.mass   = mass
        # error
        self.drift  = drift
        self.offset = offset
        self.rand   = rand
        # time
        self.t_cal  = t_cal
        self.time   = 0
        # type
        self.type   = type

    def calibrate(self):
        # reassign error rate using normal dist
        self.drift  = random.gauss(0., self.sigma_sd)
        self.offset = random.gauss(0., self.sigma_so)
        # reset time
        self.time   = 0

    # make
    def measure(self, val_in: float, dt: int = 1) -> float:
        # reassign random error using normal dist
        self.rand   = random.gauss(0., self.sigma_r)
        # calibrate if time is more than t_cal
        self.time  += dt 
        if self.time >= self.t_cal:
            self.calibrate()
        # measure
        val_out     = val_in * (1. + self.drift * self.time + self.offset + self.rand)
        # return
        return val_out

class MBA:
    # definition
    def __init__(self, IKMPs: list = []):
        # mass
        self.MUF    = 0.
        # measurement
        self.IKMPs  = IKMPs
        # time
        self.time   = 0

    # reset
    def reset(self):
        self.MUF    = 0.
        self.time   = 0

    # measurement
    def measure(self, dt: int = 1):
        # update time
        self.time  += 1
        # measure
        muf     = 0.
        for ikmp in self.IKMPs:
            m_real  = ikmp.mass
            m_meas  = ikmp.measure(m_real, dt)
            if ikmp.type == 'flow':
                muf    += m_meas
            elif ikmp.type == 'storage':
                muf    += m_meas - m_real
        # update muf
        self.MUF   += muf

#######################################################
# simulation
def MUF_over_t(mba: MBA, N_t: int, dt: int = 1, name: str = 'graph.png'):
    # measuring
    MUFs    = []
    for t in range(0, N_t, dt):
        mba.measure(dt)
        MUFs.append(mba.MUF)

    # plot    
    '''fig, ax     = plt.subplots()
    ax.plot(range(N_t), MUFs)
    ax.set_xlabel('Waktu (hari)')
    ax.set_ylabel('MUF (kg)')
    plt.savefig(name, dpi = 300)
    # plt.show()  '''  

    # return
    return MUFs

def percent_of_detected(mba: MBA, N: int, N_t: int, dt: int = 1, sq: float = SQ_LEU):
    # counter
    detected    = 0
    # loop for N times
    for n in range(N):
        mba.reset()
        # measuring
        MUFs    = []
        for t in range(0, N_t, dt):
            mba.measure(dt)
            MUFs.append(mba.MUF)
            # cek apakah lebih dari SQ
            if MUFs[-1] >= sq:
                detected   += 1
                break
    
    # return
    ratio   = float(detected)/float(N)
    return(ratio)

def percent_of_final_detected(mba: MBA, N: int, N_t: int, dt: int = 1, sq: float = SQ_LEU):
    # counter
    detected    = 0
    # loop for N times
    for n in range(N):
        mba.reset()
        # measuring
        MUFs    = []
        for t in range(0, N_t, dt):
            mba.measure(dt)
            MUFs.append(mba.MUF)
            # cek apakah lebih dari SQ
        if MUFs[-1] >= sq:
            detected   += 1
    
    # return
    ratio   = float(detected)/float(N)
    return(ratio)

#######################################################
# test
if __name__ == '__main__':
    0