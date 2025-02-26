setenv NODEBUG yes

starver .DEV2
setenv STAR /star/u/maxwoo/star_lib/packages/.DEV2
source $STAR/setupDEV2.csh
setenv NODEBUG yes
starver TFG21g
#setup gcc 4.4.7
setup 64b
setenv STARFPE NO
