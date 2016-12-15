# Rmosek Installation Guide for Mac OS X

## Contributors

Nan Xiao, Lei Sun and Peter Carbonetto.

## Environment Setup

Install Xcode from Mac App Store. Install Xcode Command Line Tools by

```bash
xcode-select -p
```

if you haven't already done so.

## Install mosek

Download the most recent version of mosek for Mac from [1], then
extract it to `~/mosek`. Note we want the "MAC OSX 64 bit x86" version
here.

This could be done under GUI, or under terminal. For example:

```bash
wget http://download.mosek.com/stable/7.1.0.53/mosektoolsosx64x86.tar.bz2
tar -xvf mosektoolsosx64x86.tar.bz2
```

## Install mosek license

Apply for free personal academic license from [2]. Download the `.lic`
file, and copy it to the mosek folder, for example,
`~/mosek/mosek.lic`.

## Install Rmosek

It is **very important** to know that the next step for installing
Rmosek should be done under the **terminal**, instead of **RStudio**.

Open your terminal app, then:

```bash
R
```

```r
install.packages("Rmosek", type="source", INSTALL_opts="--no-multiarch", 
  configure.vars="PKG_MOSEKHOME=~/mosek/7/tools/platform/osx64x86 PKG_MOSEKLIB=mosek64", 
  repos="http://download.mosek.com/R/7")
```

Then after about 30 seconds' compilation, if successful, the result
output should be something like:

```
* DONE (Rmosek)
The downloaded source packages are in
    '/tmp/Rtmpat8Xg6/downloaded_packages'
```

Before using Rmosek, it is helpful to check that it works. 
Here is a small test example from the Rmosek package:

```r
require(Rmosek)
lo1       <- list()
lo1$sense <- "max"
lo1$c <- c(3,1,5,1)
lo1$A <- Matrix(c(3,1,2,0,
                  2,1,3,1,
                  0,2,0,3),
	       nrow=3, byrow=TRUE, sparse=TRUE)
lo1$bc <- rbind(blc = c(30,15,-Inf),
                buc = c(30,Inf,25))
lo1$bx <- rbind(blx = c(0,0,0,0),
                bux = c(Inf,10,Inf,Inf))
r      <- mosek(lo1)
```

Once you have confirmed that Rmosek is working, you could install
`REBayes` and `ashr` smoothly:

```r
install.packages('REBayes')
devtools::install_github('stephens999/ashr')
```

## Links

[1] [mosek download](https://www.mosek.com/resources/downloads)

[2] [mosek free personal academic
license](https://www.mosek.com/resources/academic-license)

[3] [Rmosek installation
flags](https://stephenslab.slack.com/archives/rtips/p1461621202000012)
