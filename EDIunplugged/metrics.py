# Created by Jon Zink  (jzink@astro.ucla.edu)
# Developed with python 3.7.1
#
# If you make use of this code, please cite:
# J. K. Zink et al. 2020a


import numpy as np
from astropy.stats import mad_std
from scipy import special

def fluxContam(params,delta_mag,delta_dist, photoAp, pxScale):
    """Flux Contamination test
    Look for transit contamination from nearby stars.
    
    Args:
        params (Required[object]): the transit parameters needed to assess
            the validity of the signal
    
        delta_mag (float): the difference in magnitudes between the target
            star and the potential contaminate in the Gaia G band 
    
        delta_dist (float): the distance between the potentially contaminating
            source and the target star in arc-seconds.
    
        photoAp (int): The number of pixels used in the aperture of the flux
            measurements. 
    
    Output:
        params : the modified transit parameters needed to assess
             the validity of the signal

    """

    fit_rp=params.tlsOut.rp_rs
    fit_b=params.impact
    deltDist=(np.sqrt(photoAp/np.pi)+1)*pxScale
    fluxRatio=10**(delta_mag/-2.5)
    fTotStar=1+fluxRatio*1/2*(1+special.erf((deltDist-delta_dist)/(2.55*np.sqrt(2))))
    params.flux_ratio=fTotStar  
    
    # if deltDist>20.4:
    #     params.fluxContaminationFP=True
    
    if (fit_rp*np.sqrt(fTotStar) + fit_b)>1.04:
        params.fluxContamFP=True
        
    elif (fit_rp)*np.sqrt(fTotStar)>0.3:
        params.fluxContamFP=True
    else:
        params.fluxContamFP=False
    
    return params
  
def outlierTran(params):
    
    """Outlier Detection
     
    Looks for outliers during the apparent transit, falsely causing the signal
    
    Args:
        params (Required[object]): the transit parameters needed to assess
            the validity of the signal
    
    Output:
        params : the modified transit parameters needed to assess
             the validity of the signal

    """
     
    # global params
    per=params.tlsOut.period
    tdur=params.tlsOut.duration
    t=params.tlsOut.folded_phase
    countTrans=params.tlsOut.per_transit_count
    mes=params.tlsOut.snr
    foldY=params.tlsOut.folded_y
    foldX=params.tlsOut.folded_phase
    
    Depth=1-np.mean(foldY[(foldX>=.5-tdur/per/2) & (foldX<=.5+tdur/per/2)])
    Depthsd=mad_std(foldY[(foldX>=.5-tdur/per/2) & (foldX<=.5+tdur/per/2)])
    noDepthsd=mad_std(foldY[(foldX<.5-tdur/per/2) | (foldX>.5+tdur/per/2)])
    
    if Depthsd>(0.4*mes-1.764)*noDepthsd:
        params.outlierTranFP=True
    else:
        params.outlierTranFP=False   
        
    return params         

def transMask(params):
    
    """Individual Transit Test
     
    Looks at the individual transits for apparent anomalies.
     
    Args:
        params (Required[object]): the transit parameters needed to assess
            the validity of the signal
    
    Output:
        params : the modified transit parameters needed to assess
             the validity of the signal

    """
    nTr=params.tlsOut.transit_count
    countTrans=params.tlsOut.per_transit_count
    calcTrans=params.tlsOut.distinct_transit_count
    sES=params.tlsOut.snr_per_transit
    snr=params.tlsOut.snr
    nTrOG=nTr
    sEStot=0
    sESReal=np.zeros(len(sES))
    for i in range((nTrOG)):
        if .6*np.median(countTrans)>=countTrans[i]:
            nTr=nTr-1
        else:    
            sEStot=sEStot+sES[i]**2
            sESReal[i]=sES[i]

            
    if nTr<params.minTransit:
        params.transMaskFP=True
    else:
        params.transMaskFP=False

    if (np.max(sESReal))>=(0.8*np.sqrt(3/params.minTransit)*np.sqrt(sEStot)):
        params.transMaskFP=True
    else:
        pass

    if np.sqrt(sEStot)<params.snrThreshold:
        params.transMaskFP=True
    else:
        pass                    
                
    return params 


def evenOdd(params):
    
    eoSig=params.tlsOut.odd_even_mismatch
        
    if eoSig>5:
        params.evenOddFP=True          
    else:
        params.evenOddFP=False
    
    return params


def unique(params):
    per=params.tlsOut.period
    tdur=params.tlsOut.duration
    falseAlarm=params.tlsOut.FAP
    thershold=1/np.float(len(params.tlsOut.periods))
    foldY=params.tlsOut.folded_y
    foldX=params.tlsOut.folded_phase
    #remove SE
    foldYnoSE=foldY[(foldX>tdur/per/2) | (foldX<1-tdur/per/2)]
    foldXnoSE=foldX[(foldX>tdur/per/2) | (foldX<1-tdur/per/2)]
    #remvoe signal
    foldYnoSE=foldYnoSE[(foldXnoSE<.5-tdur/per/2) | (foldXnoSE>.5+tdur/per/2)]
    foldXnoSE=foldXnoSE[(foldXnoSE<.5-tdur/per/2) | (foldXnoSE>.5+tdur/per/2)]

    signalTrue=1-np.mean(foldY[(foldX>=.5-tdur/per/2) & (foldX<=.5+tdur/per/2)])
    signal=np.zeros(len(foldYnoSE))
    
    for i in range(len(foldYnoSE)):
        if foldX[i]<tdur/per/2:
            signal[i]=np.median(foldYnoSE)
        elif foldX[i]>(np.max(foldXnoSE)-tdur/per/2):
            signal[i]=np.median(foldYnoSE)
        else:    
            signal[i]=np.mean(foldYnoSE[(foldXnoSE>=foldXnoSE[i]-tdur/per/2) & (foldXnoSE<=foldXnoSE[i]+tdur/per/2)])
    signal=1-signal
    falseAlarmFold=np.sqrt(2)*special.erfcinv(1/np.float(len(signal)))
    signalMax=np.max(signal)
    signalAvg=np.median(signal[signal<signalMax])
    signalSD=mad_std(signal[signal<signalMax])

    if np.isnan(falseAlarm):
        params.uniqueFP=True
    elif signalTrue<signalAvg+falseAlarmFold*signalSD:
        params.uniqueFP=True
    elif falseAlarm>thershold:
        params.uniqueFP=True 
    else:    
        params.uniqueFP=False
    return params
    

def secEclipse(params):
    per=params.tlsOut.period
    tdur=params.tlsOut.duration
    foldY=params.tlsOut.folded_y
    foldX=params.tlsOut.folded_phase
    depth=params.tlsOut.depth
    fit_b=params.impact
    
    seDepth=1-np.mean(foldY[(foldX<tdur/per/2) | (foldX>1-tdur/per/2)])
    seDepthsd=mad_std(foldY[(foldX<tdur/per/2) | (foldX>1-tdur/per/2)])   
    if 0.1*depth<seDepth:
        if fit_b>=.9:
            params.SE_found=True
            params.secEclipseFP=True
        else:
            params.SE_found=True
            params.secEclipseFP=False
    elif seDepth>3*seDepthsd:
        params.secEclipseFP=False
        params.SE_found=True
        
    else:
         params.secEclipseFP=False
         params.SE_found=False
    return params


def phaseCover(params):

    fit_P=params.tlsOut.period
    fit_t0=params.tlsOut.T0
    fit_tdur=params.tlsOut.duration
    t=params.tlsOut.folded_phase
    countTrans=params.tlsOut.per_transit_count
    cadence=fit_tdur/np.median(countTrans)
    
    newMask=np.where((t>=.5-fit_tdur/fit_P) & (t<=.5+fit_tdur/fit_P),False, True)
    newT=t[~newMask]
    newT=np.sort(newT)
    dt = newT[1:] - newT[:-1]
    
    #allowed gap
    def allow_tol(x):
        func=(16*x**4-8*x**2+2)*cadence
        return func
        
    xVal=(newT[:-1]+dt/2-.5)/(fit_tdur/fit_P)
    allow=allow_tol(xVal)
    
    if np.any(np.where(dt>=allow,True,False)):
        params.phaseCoverFP=True
    else:
        params.phaseCoverFP=False
        
    return params 

def tranDur(params):
    fit_P=params.tlsOut.period
    fit_t0=params.tlsOut.T0
    fit_tdur=params.tlsOut.duration
    
    if fit_P<.5:
        params.tdurFP=True
    elif fit_tdur/fit_P>.1:
        params.tranDurFP=True
    else:
        params.tranDurFP=False
        
    return params    
                   