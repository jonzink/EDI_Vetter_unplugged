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
        >>> params=EDI_Vetter.Go(params, delta_mag=2.7, delta_dist=1000, photoAp=25, telescope="Kepler")

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
               
          


def Go(params,delta_mag=float("Inf"),delta_dist=float("Inf"), photoAp=1, telescope=None):
    """Initialize All Vetting Metrics.

    The Go function runs all of the EDI-Vetter metrics on the transit signal.  
    
    Args:
        params (Required[object]): the transit parameters needed to assess
             the validity of the signal
        
        delta_mag (optional[float]): magnitude difference between target star
             and potential contaminate source in the Gaia G band
        
        delta_dist (optional[float]): distance between the target star and the
             potential contaminate source in arc-seconds
        
        photoAp (optional[int]): number of pixels used for the target aperture
    
        telescope (optional[string]): name of telescope data was collected from 
            (i.e. "Kepler", "K2", or "TESS")
    
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
        #############   Run EDI-Vetter   ########
        ############# Exoplanet Detection Indicator - Vetter ######

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
    params=ediMet.fluxContamination(params,delta_mag,delta_dist, photoAp,pxScale)
    params=ediMet.outlierTransit(params)
    params=ediMet.individual_transits(params)
    params=ediMet.even_odd_transit(params)
    params=ediMet.uniqueness_test(params)
    # params=ephemeris_wonder(params)
    params=ediMet.check_SE(params)
    # params=harmonic_test(params)
    # params=period_alias(params,cycle=True)
    params=ediMet.phase_coverage(params)
    params=ediMet.tdur_max(params)
    
    if params.fluxContaminationFP | params.outlierTransitFP | params.TransMaskFP | params.even_odd_transit_misfit | params.uniquenessFP | params.SeFP | params.phaseCoverFP | params.tdurFP :
        params.FalsePositive=True
    else:
        params.FalsePositive=False
    
    print("==========================================")
    print("            Vetting Report") 
    print("==========================================")   
    print("        Flux Contamination : " + str(params.fluxContaminationFP))        
    # print("         Too Many Outliers : " + str(params.outlierTransitFP))
    print("  Too Many Transits Masked : " + str(params.TransMaskFP))   
    print("Odd/Even Transit Variation : " + str(params.even_odd_transit_misfit))   
    print("      Signal is not Unique : " + str(params.uniquenessFP))
    print("   Secondary Eclipse Found : " + str(params.SeFP))
    # print(" Transit Mid-point Slipped : " + str(params.eph_slipFP))
    # print("     Strong Harmonic Found : " + str(params.harmonicFP))
    print("Low Transit Phase Coverage : " + str(params.phaseCoverFP) )
    print("Transit Duration Too Large : " + str(params.tdurFP))
    print("==========================================") 
    print("Signal is a False Positive : "+  str(params.FalsePositive))
             
        
    return params
    


           
