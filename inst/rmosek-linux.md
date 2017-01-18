## Rmosek installation on the RCC midway compute cluster

This note also applies to installation on a Linux computer.

### Contributors

Nan Xiao, Lei Sun and Peter Carbonetto.

### ssh login

```bash
ssh cnetid@midway.rcc.uchicago.edu
```

### Import modules

Use `module avail R` to check all available R versions. Assume your
package library folder is `~/Rlibs_3.2.1_gcc` on RCC:

```bash
module load gcc/4.9
module load R/3.2
export R_LIBS_USER=$HOME/Rlibs_3.2.1_gcc
```

The specific versions of the modules may change later, but the above
settings worked well in May 2016.

### Install mosek

Download and extract the most recent mosek version from [1]. Note we
want the "Linux 64 bit x86" version here:

```bash
wget http://download.mosek.com/stable/7.1.0.52/mosektoolslinux64x86.tar.bz2
tar -xvf mosektoolslinux64x86.tar.bz2
```

### upload mosek license

Apply for free personal academic license from [2], download the .lic
file to your computer, say `~/mosek.lic`.

Run the following `scp` command from **local terminal** to upload the
license file to the remote mosek folder on RCC:

```bash
scp ~/mosek.lic cnetid@midway.rcc.uchicago.edu:/home/cnetid/mosek/
```

### install Rmosek

It might be a little tricky to use all the correct flags [3] when
installing Rmosek, but basically using the following command you'll be
fine:

```bash
R
```
```r
install.packages("Rmosek", type="source", INSTALL_opts="--no-multiarch",
  configure.vars="PKG_MOSEKHOME=~/mosek/7/tools/platform/linux64x86 PKG_MOSEKLIB=mosek64",
  repos="http://download.mosek.com/R/7")
```

Then after about 30 seconds' compilation, if successful, the result
output should be like:

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

### useful links

[1] [mosek download](https://www.mosek.com/resources/downloads)

[2] [mosek free personal academic license](https://www.mosek.com/resources/academic-license)

[3] [Rmosek installation flags](https://stephenslab.slack.com/archives/rtips/p1461621202000012)
