# EDI-Vetter Unplugged
This is a program meant to identify false positive transit signals using TLS information. This program has been simplified from the full [EDI-Vetter](https://github.com/jonzink/EDI-Vetter) algorithm for easy implementation with TLS.

<a href="https://zenodo.org/badge/latestdoi/200920137"><img src="https://zenodo.org/badge/200920137.svg" alt="DOI"></a>   

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for research, development, and testing purposes. EDI-Vetter Unplugged was written in Python 3.4 

### Prerequisites

Several python packages are required to run this software. Here are a few:  TLS, Numpy, scipy, and astropy

EDI-Vetter Unplugged is meant to utilize the output provided by the [TLS](https://github.com/hippke/tls) software package. We remind users to cite both packages appropriately.  


## Running EDI-Vetter Unplugged in Python

Here we provide a quick example.

Begin by opening Python in the appropriate directory. 
```
$ python
```
Now import the necessary packages
```
>>> import EDIunplugged as EDI
>>> import transitleastsquares
```
Run the light curve file through TLS
```
>>> model = transitleastsquares(Time, Flux)
>>> tlsOut = model.power()
```
Now you can set up the EDI-Vetter Unplugged parameters object with the TLS output object. For a quick start you can enter:
```
>>> params=EDI.parameters(tlsOut)
```
For a more detailed analysis you can provide additional information about your search:
```
>>> params=EDI.parameters(tlsOut, limbDark=[0.48, 0.19], impact=0, snrThreshold=7, minTransit=3)
```
Here you can enter quadratic limb darkening values, the transit impact parameter, your desired SNR threshold, and/or the minimum number of transits considered for a valid detection. The default values have been listed in the example.

Now you can run all of the vetting metrics on the signal
```
>>> params=EDI.Go(params)
```
Once completed, EDI-Vetter unplugged will print out a vetting report:
```
 ___________ _____      _   _      _   _            
|  ___|  _  \_   _|    | | | |    | | | |           
| |__ | | | | | |______| | | | ___| |_| |_ ___ _ __ 
|  __|| | | | | |______| | | |/ _ \ __| __/ _ \ '__|
| |___| |/ / _| |_     \ \_/ /  __/ |_| ||  __/ |   
\____/|___/  \___/      \___/ \___|\__|\__\___|_|   Unplugged
   
==========================================
            Vetting Report
==========================================
        Flux Contamination : False
  Too Many Transits Masked : True
Odd/Even Transit Variation : False
      Signal is not Unique : True
   Secondary Eclipse Found : False
Low Transit Phase Coverage : False
Transit Duration Too Large : False
==========================================
Signal is a False Positive : True
```
In this case, the signal was not unique with the light curve and is likely to be a false positive. Additionally, the number of meaningful transits fell below the defined threshold.

| Output | Description |
| --- | --- |
| `fluxContaminationFP` | Was neighboring flux contamination contributing significantly? |
| `TransMaskFP` | Were the individual transits questionable?  |
| `even_odd_transit_misfit` | Does the signal deviate significantly between even and odd transits? |
| `uniquenessFP` | Does the signal appear similar to other signals within the light curve? |
| `SeFP` | Does the signal appear to have a secondary eclipse? |
| `phaseCoverFP` | Does the signal lack sufficient data to detect a meaningful transit? |
| `tdurFP` | Is the transit duration too long when compared to the period? |
| `FalsePositive` | Overall, does the signal appear to be a false positive? |


 You can access the suggested classification from EDI-Vetter Unplugged by using 
```
>>> print(params.FalsePositive)
True
>>> print(params.fluxContaminationFP)
False
```
Alternatively, you can enter information about a potential contaminating star by indicating the photometric aperture size in pixels ("photoAp"), the telescope collected from ("telescope"), the separation in arcseconds from target star and the contaminating source ("delta_dist"), and the difference in visual magnitude between the sources ("delta_mag"). Note: EDI-Vetter Unplugged is currently only applicable with "Kepler", "K2", and "TESS" telescope choices.

```
>>> params=EDI.Go(params, delta_mag=10, delta_dist=1000, photoAp=25, telescope="TESS")

```
It is important to note this is not the Full EDI-Vetter suite of vetting metric, but rather a large fraction that could be easily implemented alongside TLS. Thus, EDI-Vetter unplugged is likely to have a high completeness, but a lower detection reliability.