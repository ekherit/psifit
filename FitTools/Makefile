
# ��������� 3 ������ ����� �� �������
ifeq (,$(WORKDIR))
WORKDIR = $(shell pushd ../ 1>/dev/null && pwd )
endif

#�ӣ ��� �������� �������� ����� :) - ���������� ����������������
VERBOSEMAKE=""
#���� ������ -g ��� ������ �� ����� - ���� ������ � ����
NDEBUG=""

#���������� ��� ����� ������
#CC_LOCAL := gcc
#������ ��� ����� ������
#LD_LOCAL := gcc

#PACKAGES_LOCAL=""

# �������������� ����� ��� ����������� C/C++
COPTS  =   `root-config --cflags` -I./ 

# �������������� ����� ��� ����������� Fortran 
FOPTS  =  -fno-second-underscore 
# ���� ������� �������� � ��������

# �������������� ����� ��� ������� 
#LDOPTS = -lg2c $(shell $(ROOTSYS)/bin/root-config --libs --glibs)  -lMinuit -lNew
LDOPTS = -lg2c -Xlinker -rpath -Xlinker 'root-config --libdir'  -lMinuit -lNew
 #-shared #���� ������� ������������ ����������� ����������
#���� ���������� �������� CERNLIB, �� ����� ��������������� ���� ��������
CERNLIBRARY = "" 
#������ ���������, ���� ��� �� �������� ����������� ����� (��������� �� ����������)
# jetset74 mathlib graflib geant321 grafX11 packlib
CERNLIBS = mathlib graflib grafX11 packlib 

#���� ���������� ���������� ONLYBINARY, �� ���������� � ������ �����������
#ONLYBINARY=""
# ��� ���������� ����������� �����
#BINDIR=./
#����������������, ���� ������� ������� ������������ ����������
#LIB_SHARED = ""

# �������������� ���� (����������� ����� ���)
#LIB_LOCAL= 
# `cernlib jetset74 mathlib packlib` `root-config --libs`       

# ���������, ����� ��������� �� ����� ��������
# ����� ������� ��������� ������ ������������ ��� ������ ������, 
# ����� �� �������� bin - ��� ����������, ������� ������� ��� �������
# � ������������ ������� ���������� � ���������� example
BINARIES = 

# ������, �� ����� ������� ����� ������ ��� �������
# (��� ������ �� ����� �������� � ����������)
# � ����� ���������� ���� ���������� ��� ������
#templateTest_MODULES = test
#templateTest_LIBS    = KcTemplate



# ��������� ������ ����� �� �������
include $(WORKDIR)/KcReleaseTools/rules.mk








