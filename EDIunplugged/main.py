# Created by Jon Zink  (jzink@astro.ucla.edu)
# Developed with python 3.7.1
#
# If you make use of this code, please cite:
# J. K. Zink et al. 2020a


import EDIunplugged.metrics as ediMet
import EDIunplugged.version as version
import warnings

class parameters(object):
    """Initialize Vetting.

    The parameters object itself is just a container object. Different
    codes can perform module operations on the parmameters object.
    
    Args:

        
        tlsOut (Required[object]): output from TLS transit search
        

        limbDark (Optional[array]): array of the two quadratic limb darkening 
            parameters for the stellar host
        
        impact (Optional[float]): float representing the impact parameter of
            the transit search
    
        snrThreshold (Optional[float]): float representing SNR threshold limit
            required for detection
        
        minTransit (Optional[int]): number of transit required for detection 
    

    Example:
    
        # Working with the parameters
    
        >>> params=EDI_Vetter.paramaters(tlsOut, limbDark=[0.48, 0.19], impact=0, snrThreshold=7, minTransit=3)
        >>> params=EDI_Vetter.Go(params, deltaMag=2.7, deltaDist=1000, photoAp=25, telescope="Kepler", )

    """
    
        
    
    def __init__(self, tlsOut,limbDark=None, impact=None, snrThreshold=7, minTransit=3):
    
        super(parameters,self).__init__()
        self.tlsOut=tlsOut
        self.impact=impact
        self.limbDark=limbDark
        self.snrThreshold=snrThreshold
        self.minTransit=minTransit

        self.limbDark=limbDark

        if (limbDark is None):
            print("WARNING: assuming default limb darkening values")
            self.limbDark=[0.4804, 0.1867]
        if (impact is None):
            print("WARNING: assuming default impact parameter value")
            self.impact=0
               
def prettyPrint(params):
    print("""
 ___________ _____      _   _      _   _
|  ___|  _  \_   _|    | | | |    | | | |
| |__ | | | | | |______| | | | ___| |_| |_ ___ _ __
|  __|| | | | | |______| | | |/ _ \ __| __/ _ \ '__|
| |___| |/ / _| |_     \ \_/ /  __/ |_| ||  __/ |
\____/|___/  \___/      \___/ \___|\__|\__\___|_|   Unplugged
    """)
    print("Version "
        + version.EDI_VERSIONING
        + " ("
        + version.EDI_DATE
        + ")"
        )
    print("==========================================")
    print("            Vetting Report") 
    print("==========================================")   
    print("        Flux Contamination : " + str(params.fluxContamFP))        
    print("         Too Many Outliers : " + str(params.outlierTranFP))
    print("  Too Many Transits Masked : " + str(params.transMaskFP))   
    print("Odd/Even Transit Variation : " + str(params.evenOddFP))   
    print("      Signal is not Unique : " + str(params.uniqueFP))
    print("   Secondary Eclipse Found : " + str(params.secEclipseFP))
    # print(" Transit Mid-point Slipped : " + str(params.eph_slipFP))
    # print("     Strong Harmonic Found : " + str(params.harmonicFP))
    print("Low Transit Phase Coverage : " + str(params.phaseCoverFP) )
    print(" Transit Duration Too Long : " + str(params.tranDurFP))
    print("==========================================") 
    print("Signal is a False Positive : "+  str(params.FalsePositive))

          


def Go(params, deltaMag=float("Inf"), deltaDist=float("Inf"), photoAp=1, telescope=None, print=True):
    """Initialize All Vetting Metrics.

    The Go function runs all of the EDI-Vetter metrics on the transit signal.  
    
    Args:
        params (Required[object]): the transit parameters needed to assess
             the validity of the signal
        
        deltaMag (optional[float]): magnitude difference between target star
             and potential contaminate source in the Gaia G band
        
        deltaDist (optional[float]): distance between the target star and the
             potential contaminate source in arc-seconds
        
        photoAp (optional[int]): number of pixels used for the target aperture
    
        telescope (optional[string]): name of telescope data was collected from 
            (i.e., "Kepler", "K2", or "TESS")
    
        print (optional[boolean]): do you want EDI to print out the results
    
    Output:
        params : the modified transit parameters needed to assess
             the validity of the signal


    """
    if (telescope is None):
        print("WARNING: Telescope not provided, assuming Kepler pixel scale")
        telescope="Kepler"
        
    if telescope=="Kepler":
        pxScale=3.98
    elif telescope=="K2":
        pxScale=3.98
    elif telescope=="TESS":
        pxScale=21
    else:
        print("WARNING: Telescope not recognized, assuming Kepler pixel scale")
        pxScale=3.98       

    params.FalsePositive=False
    
#############   Run EDI-Vetter Unplugged  ########
############# Exoplanet Detection Indicator - Vetter ######
            
    params=ediMet.fluxContam(params, deltaMag, deltaDist, photoAp,pxScale)
    params=ediMet.outlierTran(params)
    params=ediMet.transMask(params)
    params=ediMet.evenOdd(params)
    params=ediMet.unique(params)
    # params=ephemeris_wonder(params)
    params=ediMet.secEclipse(params)
    # params=harmonic_test(params)
    # params=period_alias(params,cycle=True)
    params=ediMet.phaseCover(params)
    params=ediMet.tranDur(params)
    
    if params.fluxContamFP | params.outlierTranFP| params.transMaskFP | params.evenOddFP | params.uniqueFP | params.secEclipseFP | params.phaseCoverFP | params.tranDurFP :
        params.FalsePositive=True
    else:
        params.FalsePositive=False
        
    if print==True:
        prettyPrint(params)
                 
    return params
    


           
