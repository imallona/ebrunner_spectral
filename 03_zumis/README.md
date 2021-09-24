Steps:
0. Install the required soft, including R 4.1.1, `libhdf5-dev`, `libgit2-dev`, `libmagick++-dev`, satija's loomR etc
1. Run the indexing for mouse/genome/alien
2. Edit the yaml file, i.e. as `test.yaml` (this example has no BC whitelisting)
3. Uncompress the relevant GTF
3. Run, e.g. with the proper binaries within the $PATH (no conda), with `./zUMIs.sh  -y test.yaml -d /home/ubuntu/sandbox/zUMIs/ 2>&1 | tee -a logs/test_all.log`


Tips to compile R/install stuff

```
# R 4 for giulia start
# outside conda!
echo $PATH

mkdir -p  ~/soft/R
cd $_

wget https://cloud.r-project.org/src/base/R-4/R-4.1.1.tar.gz
tar xzvf R-4.1.1.tar.gz
cd R-4.1.1

sudo aptitude install pcre2-utils libpcre2-dev libpcre3-dev libpcre3 \
    libpcre2-8-0 libpcre2-32-0 libpcre2-posix0 pcre2-utils

apt-get build-dep r-base

## attempt to find libicu, at not avail - v60 is installed, and it looks for 58
# ./configure LIBnn=lib --with-recommended-packages=no \
#             LDFLAGS=-L/usr/lib/x86_64-linux-gnu \
#             --with-cairo --with-libpng --with-libtiff --with-jpeglib --enable-R-shlib

# ./configure LIBnn=lib --with-recommended-packages=no \
#             LDFLAGS=-L/usr/lib/x86_64-linux-gnu/ \
#             --with-cairo --with-libpng --with-libtiff --with-jpeglib --enable-R-shlib

# get ICU 58 source code & compile

mkdir -p ~/soft/icu
cd $_
wget https://github.com/unicode-org/icu/archive/refs/tags/release-58-3.tar.gz
tar -xzvf release-58-3.tar.gz

cd /home/ubuntu/soft/icu-release-58-3/icu4c/source
./configure
ln -s /usr/include/locale.h /usr/include/xlocale.h
make

cd ~/soft/R/R-4.1.1

./configure LIBnn=lib --with-recommended-packages=no \
            LDFLAGS=-L/home/ubuntu/soft/icu-release-58-3/icu4c/source/lib/ \
            --with-cairo --with-libpng --with-libtiff --with-jpeglib --enable-R-shlib


make
make install prefix=~/soft/R/R-4.1.1

ln -s ~/soft/R/R-4.1.1/bin/R ~/soft/R/R
```
