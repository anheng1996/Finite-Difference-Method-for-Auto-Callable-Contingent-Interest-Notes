import math
import numpy as np
import xlwt

S0=195.09
S_low=0
S_high=3*S0
barrier=S0*0.78    ##also the trigger
imax=459
jmax=300
dt=1/365
S_list=np.arange(S_low,S_high+0.0000001,S0/100)

def CN(r,div,sigma):
    grid=np.ones((jmax+1,imax+1))*0####with trigger
    for i in range(0,jmax+1):           ####maturity
        if S_list[i]>=S0:
            grid[i,imax]=(1000+20.375)*math.exp(-r*3/365)
        elif S_list[i]>=barrier:
            grid[i,imax]=(1000*S_list[i]/S0+20.375)*math.exp(-r*3/365)
        else:
            grid[i,imax]=(1000*S_list[i]/S0)*math.exp(-r*(3/365))
        
    list_one=np.arange(0,jmax+0.1,1)
    a_list=[0.25*(sigma**2*a**2-(r-div)*a) for a in list_one]
    a_list[0]=a_list[jmax]=0
    
    b_list=[-sigma**2*b**2/2-r/2-1/dt for b in list_one]
    b_list[0]=1
    b_list[jmax]=1
    
    c_list=[0.25*(sigma**2*c**2+(r-div)*c) for c in list_one]
    c_list[0]=c_list[jmax]=0
    
    alpha=[b_list[0]]*(jmax+1)
    for h in range (1,jmax+1):
        alpha[h]=b_list[h]-(a_list[h]*c_list[h-1])/alpha[h-1]
        
    for i in range(imax-1,-1,-1):

        d=[0]*(jmax+1)
        for m in range(1,jmax):
            d[m]=-a_list[m]*grid[m-1,i+1]-(-sigma**2*m**2/2-r/2+1/dt)*grid[m,i+1]-c_list[m]*grid[m+1,i+1]
        
        if i<=92:
            d[jmax]=(1000+20.375)*math.exp(-r*(92-i)/365)
        elif i<=186:
            d[jmax]=(1000+20.375)*math.exp(-r*(186-i)/365)
        elif i<=277:
            d[jmax]=(1000+20.375)*math.exp(-r*(277-i)/365)
        elif i<=368:
            d[jmax]=(1000+20.375)*math.exp(-r*(368-i)/365)
        else:
            d[jmax]=(1000+20.375)*math.exp(-r*(459-i)/365)
        d[0]=0
        
        S1_list=[d[0]]*(jmax+1)
        for u in range(1,jmax+1):
            S1_list[u]=d[u]-S1_list[u-1]*a_list[u]/alpha[u-1]
            
        grid[jmax,i]=S1_list[jmax]/alpha[jmax]
        for q in range(jmax-1,-1,-1):
            grid[q,i]=(S1_list[q]-c_list[q]*grid[q+1,i])/alpha[q]
        
        if (i==92):
            for j in range(0,jmax+1):
                if S_list[j]>=barrier:
                    grid[j,i]=grid[j,i]+20.375*math.exp(-r*3/365)
                else:
                    continue
        
        if (i==186) or (i==277) or (i==368):
            for j in range(0,jmax+1):
                if S_list[j]>=S0:
                    grid[j,i]=(1000+20.375)*math.exp(-r*3/365)
                elif S_list[j]>=barrier:
                    grid[j,i]=grid[j,i]+20.375*math.exp(-r*3/365)
                else:
                    continue
                    
        else:
            continue
    
    grid1=np.ones((jmax+1,imax+1))*0   ##without trigger
    for i in range(78,jmax+1):           ####maturity
        grid1[i,imax]=(1000+20.375)*math.exp(-r*3/365)
    
    for i in range(imax-1,-1,-1):
        d=[0]*(jmax+1)
        for m in range(78,jmax):
            d[m]=-a_list[m]*grid1[m-1,i+1]-(-sigma**2*m**2/2-r/2+1/dt)*grid1[m,i+1]-c_list[m]*grid1[m+1,i+1]
        
        if i<=92:
            d[jmax]=(1000+20.375)*math.exp(-r*(92-i)/365)
        elif i<=186:
            d[jmax]=(1000+20.375)*math.exp(-r*(186-i)/365)
        elif i<=277:
            d[jmax]=(1000+20.375)*math.exp(-r*(277-i)/365)
        elif i<=368:
            d[jmax]=(1000+20.375)*math.exp(-r*(368-i)/365)
        else:
            d[jmax]=(1000+20.375)*math.exp(-r*(459-i)/365)
        d[78]=grid[78,i]
        
        S1_list=[d[78]]*(jmax+1)
        for u in range(78,jmax+1):
            S1_list[u]=d[u]-S1_list[u-1]*a_list[u]/alpha[u-1]
            
        grid1[jmax,i]=S1_list[jmax]/alpha[jmax]
        for q in range(jmax-1,-1,-1):
            grid1[q,i]=(S1_list[q]-c_list[q]*grid1[q+1,i])/alpha[q]
        
        if (i==92):
            for j in range(0,jmax+1):
                if S_list[j]>=barrier:
                    grid1[j,i]=grid1[j,i]+20.375*math.exp(-r*3/365)
                else:
                    continue
        
        if (i==186) or (i==277) or (i==368):
            for j in range(0,jmax+1):
                if S_list[j]>=S0:
                    grid1[j,i]=(1000+20.375)*math.exp(-r*3/365)
                elif S_list[j]>=barrier:
                    grid1[j,i]=grid1[j,i]+20.375*math.exp(-r*3/365)
                else:
                    continue
                    
        else:
            continue
    
    return(grid[100,0])

value=CN(0.025756,0.01644,0.24686)