#ifndef ENUMS_H
#define ENUMS_H

enum CentType
{
   Cent0to5 = 0, Cent5to10 = 1, Cent10to20 = 2, Cent20to30 = 3,
   Cent30to40 = 4, Cent40to50 = 5, Cent50to70 = 6, Cent70to100 = 7
};
enum AbsEtaType
{
   AbsEta0p0to1p0 = 0, AbsEta1p0to2p0 = 1
};
enum RType
{
   R0p3 = 0, R0p4 = 1, R0p6 = 2, R0p8 = 3, R1p0 = 4
};
enum FlowType
{
   NoFlow = 0,
   FlowDefaultInRho = 1
};


#endif
