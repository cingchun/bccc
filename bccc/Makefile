# create ta
obts = pub.o \
       ta/ta_.o \
       ta/ta_A10B11.o \
       ta/ta_A10B12.o \
       ta/ta_A11B10.o \
       ta/ta_A11B9.o \
       ta/ta_A12B10.o \
       ta/ta_A12B9.o \
       ta/ta_A13B7.o \
       ta/ta_A13B8.o \
       ta/ta_A14B7.o \
       ta/ta_A14B8.o \
       ta/ta_A15B6.o \
       ta/ta_A1B1.o \
       ta/ta_A1B2.o \
       ta/ta_A1B3.o \
       ta/ta_A2B1.o \
       ta/ta_A2B2.o \
       ta/ta_A2B3.o \
       ta/ta_A3B1.o \
       ta/ta_A3B2.o \
       ta/ta_A3B3.o \
       ta/ta_A4B5.o \
       ta/ta_A5B4.o \
       ta/ta_A6B15.o \
       ta/ta_A7B13.o \
       ta/ta_A7B14.o \
       ta/ta_A8B13.o \
       ta/ta_A8B14.o \
       ta/ta_A9B11.o \
       ta/ta_A9B12.o


ta.so : $(obts)
    cc -std=c++11 -qopenmp -fPIC -shared -o ta.so $(obts)

$(obts) : pub.h ta/ta.h
    cc -o