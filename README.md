# EDI-Vetter unplugged
This is a program meant identify false positive transit signal using TLS info. This program has been simplified from the full EDI-Vetter algorithm for easy implementation with TLS.

<a href="https://zenodo.org/badge/latestdoi/200920137"><img src="https://zenodo.org/badge/200920137.svg" alt="DOI"></a>   

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for research, development, and testing purposes. EDI-Vetter unnplugged was written in Python 3.4 

### Prerequisites

Several python packages are required to run this software. Here are a few:  Numpy,  scipy, and astropy

EDI-Vetter unplugged is meant to utilize the output dataframe provided by the [TLS](https://github.com/hippke/tls) software package. We have included a copy in this repositories, but remind users to cite appropriately.  




## Running EDI-Vetter unplugged in Python

Here we provide a quick example using the light curve of K2-138.

Begin by opening Python in the appropriate directory. 
```
$ python
```
Now import the necessary packages
```
>>> import EDI_Vetter_unplugged
>>> import transitleastsquares
```
Run the light curve file through TLS
```
>>> model = transitleastsquares(Time, Flux)
>>> tlsOut = model.power()
```
Now you can set up the EDI-Vetter unplugged parameters object with the TLS output dataframe 
```
>>> params=EDI_Vetter_unplugged.parameters(tlsOut, limbDark=[0.48, 0.19], impact=0, snrThreshold=7, minTransit=3)
```
Here you can enter quadratic limb darkening values, the transit impact parameter, your desired SNR threshold, and/or the minimum number of transits to be considered valid. The default values have been listed in the example.
Now you can run all of the vetting metrics on the signal
```
>>> params=EDI_Vetter_unplugged.Go(params,delta_mag=10,delta_dist=1000, photoAp=41)
```
You can enter a potential contaminating star by indicating the photometric aperture size in pixels ("photoAp"), the separation in arcseconds ("delta_dist"), and the visual magnitude of the contaminating source ("delta_mag"). Note: Kepler pixels are assumed here. You do not need to enter any of these parameters if you are not concerned about flux contamination.