# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 19:25:47 2023

@author: sun'yi'xin
"""
import gamspy.math as gams_math
from gamspy import Card
from gamspy import Container
from gamspy import Equation
from gamspy import Model
from gamspy import Ord
from gamspy import Parameter
from gamspy import Set
from gamspy import Sum
from gamspy import Variable


def main():
    m = Container(delayed_execution=True)

    # SETS #
    t = Set(
        m,
        name="t",
        records=[str(t) for t in range(1, 61)],
        description="time periods",
    )
    tfirst = Set(
        m, name="tfirst", domain=[t], description="first interval (t0)"
    )
    tlast = Set(m, name="tlast", domain=[t], description="last intervat [T]")
    tnotlast = Set(
        m, name="tnotlast", domain=[t], description="all intervals but last"
    )
    lag10 =  Set(m, name="lag10", domain=[t], description="last 10 years")
    
    lag10[t].where[Ord(t) >= Card(t)-10] = True
    tfirst[t].where[Ord(t) == 1] = True
    tlast[t].where[Ord(t) == Card(t)] = True
    tnotlast[t] = ~tlast[t]
    

    # SCALARS #
    
    tstep = Parameter(m, name="tstep", records=5, description="Years per Period")
    
    #If optimal control
    ifopt = Parameter(m, name="ifopt", records=0, description="Indicator where optimized is 1 and base is 0")
    
    #Preferences
    elasmu = Parameter(m, name="elasmu", records=1.45, description="Elasticity of marginal utility of consumption")
    prstp = Parameter(m, name="prstp", records=0.015, description="Initial rate of social time preference per year")
    
    #Population and technology
    gama = Parameter(m, name="gama", records=0.3, description="Capital elasticity in production function")
    pop0 = Parameter(m, name="pop0", records=6838, description="Initial world population (millions)")
    popadj = Parameter(m, name="popadj", records=0.134, description="Growth rate to calibrate to 2050 pop projection")
    popasym = Parameter(m, name="popasym", records=10500, description="Asymptotic population (millions)")
    dk = Parameter(m, name="dk", records=0.1, description="Depreciation rate on capital (per year)")
    q0 = Parameter(m, name="q0", records=63.69, description="Initial world gross output (trill 2005 USD)")
    k0 = Parameter(m, name="k0", records=135, description="Initial capital value (trill 2005 USD)")
    a0 = Parameter(m, name="a0", records=3.80, description="Initial level of total factor productivity")
    ga0 = Parameter(m, name="ga0", records=0.079, description="Initial growth rate for TFP per 5 years")
    dela = Parameter(m, name="dela", records=0.006, description="Decline rate of TFP per year")

    #Emissions parameters
    gsigma1 = Parameter(m, name="gsigma1", records=-0.01, description="Initial growth of sigma (per year)")
    dsig = Parameter(m, name="dsig", records=-0.001, description="Decline rate of decarbonization (per period)")
    eland0 = Parameter(m, name="eland0", records=3.3, description="Carbon emissions from land 2010 (GtCO2 per year)")
    deland = Parameter(m, name="deland", records=0.2, description="Decline rate of land emissions (per period)")
    e0 = Parameter(m, name="e0", records=33.61, description="Industrial emissions 2010 (GtCO2 per year)")
    miu0 = Parameter(m, name="miu0", records=.039, description="Initial emissions control rate for base case 2010")
    
    ##Carbon cycle
    #Initial Conditions
    mat0 = Parameter(m, name="mat0", records=830.4, description="Initial Concentration in atmosphere 2010 (GtC)")
    mu0 = Parameter(m, name="mu0", records=1527., description="Initial Concentration in upper strata 2010 (GtC)")
    ml0 = Parameter(m, name="ml0", records=10010., description="Initial Concentration in lower strata 2010 (GtC)")
    MATEQ = Parameter(m, name="MATEQ", records=588, description="Equilibrium concentration atmosphere  (GtC)")
    mueq = Parameter(m, name="mueq", records=1350, description="Equilibrium concentration in upper strata (GtC)")
    mleq = Parameter(m, name="mleq", records=10000, description="Equilibrium concentration in lower strata (GtC)")
    
    #Flow paramaters
    b12 = Parameter(
        m, name="b12", 
        records=.088, 
        description="Carbon cycle transition matrix"
    )
    b23 = Parameter(
        m, name="b23", 
        records=0.00250, 
        description="Carbon cycle transition matrix"
    )

    #These are for declaration and are defined later
    b11 = Parameter(
        m, name="b11", 
        description="Carbon cycle transition matrix"
    )
    b21 = Parameter(
        m, name="b21", 
        description="Carbon cycle transition matrix"
    )
    b22 = Parameter(
        m, name="b22", 
        description="Carbon cycle transition matrix"
    )
    b32 = Parameter(
        m, name="b32", 
        description="Carbon cycle transition matrix"
    )
    b33 = Parameter(
        m, name="b33", 
        description="Carbon cycle transition matrix"
    )
    sig0 = Parameter(
        m, name="sig0", 
        description="Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)"
    )
     
    #Climate model parameters
    t2xco2 = Parameter(
        m, name="t2xco2", records=2.9, description="Equilibrium temp impact (oC per doubling CO2)"
    )
    fex0 = Parameter(
        m, name="fex0", records=0.25, description="2010 forcings of non-CO2 GHG (Wm-2)"
    )
    fex1 = Parameter(
        m, name="fex1", records=0.70, description="2100 forcings of non-CO2 GHG (Wm-2)"
    )
    tocean0 = Parameter(
        m, name="tocean0", records=.0068, description="Initial lower stratum temp change (C from 1900)"
    )
    tatm0 = Parameter(
        m, name="tatm0", records=0.80, description="Initial atmospheric temp change (C from 1900)"
    )
    
    c10 = Parameter(
        m, name="c10", records=0.098, description="Initial climate equation coefficient for upper level"
    )
    c1beta = Parameter(
        m, name="c1beta", records=0.01243, description="Regression slope coefficient(SoA~Equil TSC)"
    )
    
    c1 = Parameter(
        m, name="c1", records=0.098, description="Climate equation coefficient for upper level"
    )
    c3 = Parameter(
        m, name="c3", records=0.088, description="Transfer coefficient upper to lower stratum"
    )
    c4 = Parameter(
        m, name="c4", records=0.025, description="Transfer coefficient for lower level"
    )
    fco22x = Parameter(
        m, name="fco22x", records=3.8, description="Forcings of equilibrium CO2 doubling (Wm-2)"
    )
    
    #Climate damage parameters
    a10 = Parameter(
        m, name="a10", records=0, description="Initial damage intercept"
    )
    a20 = Parameter(
        m, name="a20", records=0.00267, description="Initial damage quadratic term"
    )
    a1 = Parameter(
        m, name="a1", records=0, description="Damage intercept "
    )
    a2 = Parameter(
        m, name="a2", records=0.00267, description="Damage quadratic term"
    )
    a3 = Parameter(
        m, name="a3", records=2.00, description="Damage exponent"
    )
    
    #Abatement cost
    expcost2 = Parameter(
        m, name="expcost2", records=2.8, description="Exponent of control cost function"
    )
    pback = Parameter(
        m, name="pback", records=344, description="Cost of backstop 2005$ per tCO2 2010"
    )
    gback = Parameter(
        m, name="gback", records=.025, description="Initial cost decline backstop cost per period"
    )
    limmiu = Parameter(
        m, name="limmiu", records=1.2, description="Upper limit on control rate after 2150"
    )
    tnopol = Parameter(
        m, name="tnopol", records=45, description="Period before which no emissions controls base"
    )
    cprice0 = Parameter(
        m, name="cprice0", records=1.0, description="Initial base carbon price (2005$ per tCO2)"
    )
    gcprice = Parameter(
        m, name="gcprice", records=.02, description="Growth rate of base carbon price per year"
    )
    
    #Participation parameters
    periodfullpart = Parameter(
        m, name="periodfullpart", records=21, description="Period at which have full participation"
    )
    partfract2010 = Parameter(
        m, name="partfract2010", records=1, description="Fraction of emissions under control in 2010"
    )
    partfractfull = Parameter(
        m, name="partfractfull", records=1, description="Fraction of emissions under control at full time"
    )
    
    #Availability of fossil fuels
    fosslim = Parameter(
        m, name="fosslim", records=6000, description=" Maximum cumulative extraction fossil fuels (GtC)"
    )

    ##Scaling and inessential parameters
    #Note that these are unnecessary for the calculations but are for convenience
    scale1 = Parameter(
        m, name="scale1", records=0.016408662, description="Multiplicative scaling coefficient"
    )
    scale2 = Parameter(
        m, name="scale2", records=-3855.106895, description="Additive scaling coefficient"
    )
    
    # PARAMETERS #
    l = Parameter(
        m, name="l", domain=[t], description="labor (production input)"
    )
    al = Parameter(
        m, name="al", domain=[t], description="Level of total factor productivity"
    )
    sigma = Parameter(
        m, name="sigma", domain=[t], description="CO2-equivalent-emissions output rati"
    )
    rr = Parameter(
        m, name="rr", domain=[t], description="Average utility social discount rate"
    )
    ga = Parameter(
        m, name="ga", domain=[t], description="Growth rate of productivity from"
    )
    forcoth = Parameter(
        m, name="forcoth", domain=[t], description="Exogenous forcing for other greenhouse gases"
    )
    gl = Parameter(
        m, name="gl", domain=[t], description="Growth rate of labor"
    )
    gcost1 = Parameter(
        m, name="gcost1", domain=[t], description="Growth of cost factor"
    )
    gsig = Parameter(
        m, name="gsig", domain=[t], description="Change in sigma (cumulative improvement of energy efficiency)"
    )
    etree = Parameter(
        m, name="etree", domain=[t], description="Emissions from deforestation"
    )
    cost1 = Parameter(
        m, name="cost1", domain=[t], description="Adjusted cost for backstop"
    )
    partfract = Parameter(
        m, name="partfract", domain=[t], description="Fraction of emissions in control regime"
    )
    lam = Parameter(
        m, name="lam", description="Climate model parameter"
    )
    gfacpop = Parameter(
        m, name="gfacpop", domain=[t], description="Growth factor population"
    )
    pbacktime = Parameter(
        m, name="pbacktime", domain=[t], description="Backstop price"
    )
    optlrsav = Parameter(
        m, name="optlrsav", description="Optimal long-run savings rate used for transversality"
    )
    scc = Parameter(
        m, name="scc", domain=[t], description="Social cost of carbon"
    )
    cpricebase = Parameter(
        m, name="cpricebase", domain=[t], description="Carbon price in base case"
    )
    photel = Parameter(
        m, name="photel", domain=[t], description="Carbon Price under no damages (Hotelling rent condition)"
    )
    
    #Parameters for long-run consistency of carbon cycle
    b11 = 1 - b12
    b21 = b12 * MATEQ / mueq
    b22 = 1 - b21 - b23
    b32 = b23 * mueq / mleq
    b33 = 1 - b32
  
    #Further definitions of parameters
    sig0 = e0 / ( q0 * (1-miu0) )
    lam = fco22x/ t2xco2
    l[tfirst] = pop0 
    for i in tnotlast:
        l[t.lead(1)] = l[t]
        l[t.lead(1)]=l[t]*(popasym/(l[t]+0.0000001))**popadj

    ga[t]=ga0*gams_math.exp(-dela*5*((t.val-1)))
    al[tfirst] = a0
    for i in tnotlast:
        al[t.lead(1)]=al[t]/((1-ga[t]))
    
    gsig[tfirst]=gsigma1
    for i in tnotlast:
        gsig[t.lead(1)]=gsig[t]*((1+dsig)**tstep)

    sigma[tfirst]=sig0
    for i in tnotlast:
        sigma[t.lead(1)]=sigma[t]*gams_math.exp(gsig[t]*tstep)

    pbacktime[t]= pback * (1-gback) ** ((t.val-1))
    cost1[t] = pbacktime[t]*sigma[t]/expcost2/1000
    etree[t] = eland0*(1-deland)**((t.val-1))
    rr[t] = 1/((1+prstp)**(tstep*(t.val-1)))
    if  Ord(t) <= 18:
        forcoth[t] = fex0 + (1/18)*(fex1-fex0)*(t.val-1)
    else:
        forcoth[t] = fex0 + (fex1-fex0)
    optlrsav = (dk + .004)/(dk + .004*elasmu + prstp)*gama
    if Ord(t) > periodfullpart:
        partfract[t] = partfractfull
    if Ord(t) < periodfullpart+1:
        partfract[t] = partfract2010+(partfractfull-partfract2010)*(Ord(t)-1)/periodfullpart

    partfract[tfirst]= partfract2010

    #Transient TSC Correction ("Speed of Adjustment Parameter")
    c1 =  c10 + c1beta*(t2xco2-2.9)

    #Base Case      Carbon Price
    cpricebase[t]= cprice0*(1+gcprice)**(5*(t.val-1))
    
    
    
    

    # VARIABLES #
    MIU = Variable(m, name="MIU", domain=[t], description="Emission control rate GHGs")
    FORC = Variable(m, name="FORC", domain=[t], description="Increase in radiative forcing (watts per m2 from 1900)")
    FORC1 = Variable(m, name="FORC1", domain=[t], description="Increase in radiative forcing adjust")
    TATM = Variable(m, name="TATM", domain=[t], description="Increase temperature of atmosphere (degrees C from 1900)")
    TOCEAN = Variable(m, name="TOCEAN", domain=[t], description="Increase temperatureof lower oceans (degrees C from 1900)")
    MAT = Variable(m, name="MAT", domain=[t], description="Carbon concentration increase in atmosphere (GtC from 1750)")
    MU = Variable(m, name="MU", domain=[t], description="Carbon concentration increase in shallow oceans (GtC from 1750)")
    ML = Variable(m, name="ML", domain=[t], description="Carbon concentration increase in lower oceans")
    E = Variable(m, name="E", domain=[t], description="Total CO2 emissions")
    EIND = Variable(m, name="EIND", domain=[t], description="Industrial emissions")
    C = Variable(m, name="C", domain=[t], description="Consumption")
    K = Variable(m, name="K", domain=[t], description="Capital stock")
    CPC = Variable(m, name="CPC", domain=[t], description="Per capita consumption")
    I = Variable(m, name="I", domain=[t], description="Investment")
    S = Variable(m, name="S", domain=[t], description="Gross savings rate")
    RI = Variable(m, name="RI", domain=[t], description="Real interest rate")
    Y = Variable(m, name="Y", domain=[t], description="Gross world product net of abatement and damages")
    YGROSS = Variable(m, name="YGROSS", domain=[t], description="Gross world product GROSS of abatement and damages")
    YNET = Variable(m, name="YNET", domain=[t], description="Output net of damages equation")
    DAMAGES = Variable(m, name="DAMAGES", domain=[t], description="Damages")
    DAMFRAC = Variable(m, name="DAMFRAC", domain=[t], description="Damages as fraction of gross output")
    ABATECOST = Variable(m, name="ABATECOST", domain=[t], description="Cost of emissions reductions")
    MCABATE = Variable(m, name="MCABATE", domain=[t], description="Marginal cost of abatement")
    CCA = Variable(m, name="CCA", domain=[t], description="Cumulative industrial carbon emissions")
    PERIODU = Variable(m, name="PERIODU", domain=[t], description="One period utility function")
    CPRICE = Variable(m, name="CPRICE", domain=[t], description="Carbon price")
    CEMUTOTPER = Variable(m, name="CEMUTOTPER", domain=[t], description="Period utility")
    UTILITY = Variable(m, name="UTILITY", description="Welfare function")

    # EQUATIONS #
    #Emissions and Damages
    EEQ = Equation(
        m, name="EEQ", type="regular", domain=[t], description="Emissions equation"
    )
    EINDEQ = Equation(
        m,
        name="EINDEQ",
        type="regular",
        domain=[t],
        description="Industrial emissions",
    )
    CCACCA = Equation(
        m,
        name="CCACCA",
        type="regular",
        domain=[t],
        description="Cumulative carbon emissions",
    )
    FORCE = Equation(
        m,
        name="FORCE",
        type="regular",
        domain=[t],
        description="Radiative forcing equation",
    )
    FORCE1 = Equation(
        m,
        name="FORCE1",
        type="regular",
        domain=[t],
        description="Radiative forcing equation1",
    )
    DAMFRACEQ = Equation(
        m,
        name="DAMFRACEQ",
        type="regular",
        domain=[t],
        description="Equation for damage fraction",
    )
    DAMEQ = Equation(
        m,
        name="DAMEQ",
        type="regular",
        domain=[t],
        description="Damage equation",
    )
    ABATEEQ = Equation(
        m,
        name="ABATEEQ",
        type="regular",
        domain=[t],
        description="Cost of emissions reductions equation",
    )
    MCABATEEQ = Equation(
        m,
        name="MCABATEEQ",
        type="regular",
        domain=[t],
        description="Equation for MC abatement",
    )
    CARBPRICEEQ = Equation(
        m,
        name="CARBPRICEEQ",
        type="regular",
        domain=[t],
        description="Carbon price equation from abatement",
    )
    #Climate and carbon cycle
    MMAT = Equation(
        m,
        name="MMAT",
        type="regular",
        domain=[t],
        description="Atmospheric concentration equation",
    )
    MMU = Equation(
        m,
        name="MMU",
        type="regular",
        domain=[t],
        description="Shallow ocean concentration",
    )
    MML = Equation(
        m,
        name="MML",
        type="regular",
        domain=[t],
        description="Lower ocean concentration",
    )
    TATMEQ = Equation(
        m,
        name="TATMEQ",
        type="regular",
        domain=[t],
        description="Temperature-climate equation for atmosphere",
    )
    TOCEANEQ = Equation(
        m,
        name="TOCEANEQ",
        type="regular",
        domain=[t],
        description="Temperature-climate equation for lower oceans",
    )
    #Economic variables
    YGROSSEQ = Equation(
        m,
        name="YGROSSEQ",
        type="regular",
        domain=[t],
        description="Output gross equation",
    )
    YNETEQ = Equation(
        m,
        name="YNETEQ",
        type="regular",
        domain=[t],
        description="Output net of damages equation",
    )
    YY = Equation(
        m,
        name="YY",
        type="regular",
        domain=[t],
        description="Output net equation",
    )
    CC = Equation(
        m,
        name="CC",
        type="regular",
        domain=[t],
        description="Consumption equation",
    )
    CPCE = Equation(
        m,
        name="CPCE",
        type="regular",
        domain=[t],
        description="Per capita consumption definition",
    )
    SEQ = Equation(
        m,
        name="SEQ",
        type="regular",
        domain=[t],
        description="Savings rate equation",
    )
    KK = Equation(
        m,
        name="KK",
        type="regular",
        domain=[t],
        description="Capital balance equation",
    )
    RIEQ = Equation(
        m,
        name="RIEQ",
        type="regular",
        domain=[t],
        description="Interest rate equation",
    )
    #Utility
    CEMUTOTPEREQ = Equation(
        m,
        name="CEMUTOTPEREQ",
        type="regular",
        domain=[t],
        description="Period utility",
    )
    PERIODUEQ = Equation(
        m,
        name="PERIODUEQ",
        type="regular",
        domain=[t],
        description="Instantaneous utility function equation",
    )
    UTIL = Equation(
        m,
        name="UTIL",
        type="regular",
        description="Objective function",
    )


    ##Equations of the model
    #Emissions and Damages
    EEQ[t] = E[t] == EIND[t] + etree[t]
    EINDEQ[t]= EIND[t] == sigma[t] * YGROSS[t] * (1-(MIU[t]))
    CCACCA[tnotlast[t]] = CCA[t.lead(1)] == CCA[t]+ EIND[t]*5/3.666
    FORCE1[t] = FORC1[t] == gams_math.log(MAT[t])
    FORCE[t] = FORC[t]  == fco22x * FORC1[t]/588.000/gams_math.log(2) + forcoth[t.val]
    DAMFRACEQ[t] = DAMFRAC[t] == (a1*TATM[t])+(a2*TATM[t]**a3) 
    DAMEQ[t] = DAMAGES[t] == YGROSS[t] * DAMFRAC[t]
    ABATEEQ[t] = ABATECOST[t]  == YGROSS[t] * cost1[t] * (MIU[t]**expcost2)* (partfract[t]**(1-expcost2))
    MCABATEEQ[t] = MCABATE[t]  == pbacktime[t] * MIU[t]**(expcost2-1)
    CARBPRICEEQ[t] = CPRICE[t]  == pbacktime[t] * (MIU[t])**(expcost2-1)/partfract[t]**(expcost2-1)

    #Climate and carbon cycle
    MMAT[tnotlast[t]] = MAT[t.lead(1)]  == MAT[t]*b11 + MU[t]*b21 + (E[t]*(5/3.666))
    MML[tnotlast[t]] =  ML[t.lead(1)]   == ML[t]*b33  + MU[t]*b23
    MMU[tnotlast[t]]  = MU[t.lead(1)]   == MAT[t]*b12 + MU[t]*b22 + ML[t]*b32
    TATMEQ[tnotlast[t]] = TATM[t.lead(1)]  == TATM[t] + c1 * ((FORC[t.lead(1)]-(fco22x/t2xco2)*TATM[t])-(c3*(TATM[t]-TOCEAN[t])))
    TOCEANEQ[tnotlast[t]] = TOCEAN[t.lead(1)]  == TOCEAN[t] + c4*(TATM[t]-TOCEAN[t])

    #Economic variables
    YGROSSEQ[t] =  YGROSS[t]  == (al[t]*(l[t]/1000)**(1-gama))*(K[t]**gama)
    YNETEQ[t] = YNET[t]   == YGROSS[t]*(1-DAMFRAC[t])
    YY[t]  =  Y[t]  == YNET[t] - ABATECOST[t]
    CC[t]  =  C[t]  == Y[t] - I[t]
    CPCE[t] = CPC[t]  == 1000 * C[t] / (l[t]+0.0000001)
    SEQ[t] = I[t] == S[t] * Y[t]
    KK[tnotlast[t]] = K[t.lead(1)]  == (1-dk)**tstep * K[t] + tstep * I[t]
    RIEQ[tnotlast[t]]  = RI[t]  == (1+prstp) * (CPC[t.lead(1)]/(CPC[t]+0.00000001))**(elasmu/tstep) - 1

    #Utility
    CEMUTOTPEREQ[t] =  CEMUTOTPER[t]  == PERIODU[t] * l[t] * rr[t]
    PERIODUEQ[t] = PERIODU[t] == ((C[t]*1000/(l[t]+0.00000001))**(1-elasmu)-1)/(1-elasmu)-1
    UTIL[...] =  UTILITY  == tstep * scale1 * Sum(t,  CEMUTOTPER[t]) + scale2

    #Resource limit
    CCA.up[t]  = fosslim
    #Control rate limits
    MIU.up[t] = limmiu*partfract[t]
    MIU.up[t] = 1.2

    #Upper and lower bounds for stability
    K.lo[t]         = 1
    MAT.lo[t]       = 10
    MU.lo[t]        = 100
    MIU.lo[t]       = 0.000001
    ML.lo[t]        = 1000
    C.lo[t]         = 2
    TOCEAN.up[t]    = 20
    TOCEAN.lo[t]    = -1
    TATM.up[t]      = 40
    CPC.lo[t]       = .01

    #control variables
    #Set savings rate for steady state for last 10 periods
    S.fx[lag10[t]] = optlrsav

    #Initial conditions
    CCA.fx[tfirst]    = 90
    K.fx[tfirst]      = k0
    MAT.fx[tfirst]    = mat0
    MU.fx[tfirst]     = mu0
    #MIU.fx[tfirst]    = miu0
    ML.fx[tfirst]     = ml0
    TATM.fx[tfirst]   = tatm0
    TOCEAN.fx[tfirst]  = tocean0
    
    co2 = Model(
        m,
        name="co2",
        equations=m.getEquations(),
        problem="nlp",
        sense="MAX",
        objective=UTILITY,
    )
    
    if ifopt == 0:
       a2 = 0
       co2.solve()


    if ifopt == 1:
        a2 = a20
        MIU.fx[tfirst] = miu0
        co2.solve()


    print("Objective Function Value:  ", round(UTILITY.toValue(), 4))

    # Solution visualization
    # ----------------------

    rep = Parameter(m, name="rep", domain=[t, "*"])

    rep[t, "C[t]"] = C.l[t]
    rep[t, "K[t]"] = K.l[t]
    rep[t, "EIND[t]"] = EIND.l[t]
    rep[t, "MAT[t]"] = MAT.l[t]
    rep[t, "TATM[t]"] = TATM.l[t]
    rep[t, "Y[t]"] = Y.l[t]
    rep[t, "CPC[t]"] = CPC.l[t]
    rep[t, "MIU[t]"] = MIU.l[t]
    rep[t, "TATM[t]"] = TATM.l[t]

    print("Solution:\n", rep.pivot().round(3))
    # End CO2

if __name__ == "__main__":
    main()
