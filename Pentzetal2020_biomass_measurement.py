import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


X=[]
Y=[]
numbins=10
rep=[]
rep_Vol_weighted=[]
strain_list=[]
rep_vol_weighted=[]


data1 = np.genfromtxt("3AUG16_analysis_good.csv", dtype=float, delimiter=',')
data2 = np.genfromtxt("3AUG16_analysis_good.csv", dtype=str, delimiter=',')
Vol=data1[1:,5]
Volsnow=data1[1:,7]
Volfloc=data1[1:,6]
Vol=np.array(Vol)
Volsnow=np.array(Volsnow)
Volfloc=np.array(Volfloc)
LogVol=np.log(Vol)
LogVolsnow=np.log(Volsnow)
LogVolfloc=np.log(Volfloc)
strain=data2[1:,4]
volumes=zip(Vol, LogVol, LogVolsnow, LogVolfloc)
#putting it into a dataframe so that we can index on strain label
dataframe=pd.DataFrame(volumes, index=strain)
#print dataframe

#calculating the number of unique labels
strain=data2[1:,4]
strain_noDupes = []
[strain_noDupes.append(i) for i in strain if not strain_noDupes.count(i)]
strain=str(strain)
print("number of strains"),len(strain_noDupes)
"""#checking that indexing is working right
print dataframe.loc[strain_noDupes[14]]
tempslice=dataframe.loc[strain_noDupes[14]]
tempslice=np.array(tempslice)
print tempslice[1:,1]
#looks good!"""

colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(strain_noDupes))])

for i in range(0,len(strain_noDupes)):
    print("computing strain"), strain_noDupes[i]
    binmeans=[]
    biomass=[]
    snow_biomass=[]
    floc_biomass=[]
    binvec=[]
    biomass_vec=[]
    total_biomass=[]
    snow_percent_overall=[]
    floc_percent_overall=[]
    tempslice=dataframe.loc[strain_noDupes[i]]
    tempslice=np.array(tempslice)
    tempVol=tempslice[0:,0]
    tempLogVol=tempslice[0:,1]
    tempLogVolsnow=tempslice[0:,2]
    tempLogVolfloc=tempslice[0:,3]
    total_biomass=sum(tempslice[0:,0])
    binsize=(max(tempLogVol)-min(tempLogVol))/(numbins) #don't have to do logs but might be good for something that is long-tailed
    for d in range(0,(numbins+1)):
        binvec.append(min(tempLogVol)+binsize*d)
    for h in range(0,len(binvec)-1):
        temp_size=[]
        temp_cluster=[]
        temp_snow=[]
        temp_floc=[]
        percent_snow=[]
        percent_floc=[]
        snow_biomass_overall=[]
        floc_biomass_overall=[]
        for a in range(0,(len(tempLogVol))):
            if tempLogVol[a]>=binvec[h] and tempLogVol[a]<=binvec[h+1]:
                temp_size.append(math.exp(tempLogVol[a])*100/total_biomass) #percentage of total biomass represented by that cluster
                temp_cluster.append(math.exp(tempLogVol[a])) #appending whole cluster information
                temp_snow.append(math.exp(tempLogVolsnow[a])) #appending snowflake biomass of that cluster
                temp_floc.append(math.exp(tempLogVolfloc[a])) #appending flocbiomass biomass of that cluster
        if sum(temp_cluster) != 0: #calculating percent of that bin that is snow or floc
            percent_snow=(sum(temp_snow)*100/sum(temp_cluster))
            percent_floc=(sum(temp_floc)*100/sum(temp_cluster))
        else: #for bins that don't have any clusters in it (can't divide by 0)
            percent_snow=0
            percent_floc=0
        snow_biomass.append(percent_snow) #appending percentage snow information for bins
        floc_biomass.append(percent_floc) #appending percentage floc information for bins
        snow_biomass_overall=(sum(temp_snow)*100/total_biomass)
        floc_biomass_overall=(sum(temp_floc)*100/total_biomass)
        snow_percent_overall.append(snow_biomass_overall)
        floc_percent_overall.append(floc_biomass_overall)
        biomass.append(sum(temp_size))
    for e in range(0,len(binvec)-1):
        binmeans.append((binvec[e]+binvec[e+1])/2) 
    rep_vol_weighted.append(np.average(binmeans, weights=biomass))
    print("binmeans"), binmeans
    print("percent_snow"), snow_biomass
    print("percent_floc"), floc_biomass
    print("biomass"), biomass
    plt.plot(binmeans, snow_biomass, alpha=.9)
    plt.plot(binmeans, floc_biomass, alpha=.9)
    plt.xlabel('Cluster size (um^2)')
    plt.ylabel('% biomass')  
    plt.xlim([0,14])
    plt.tight_layout()
    plt.savefig("%s snow_floc_biomass.pdf" % strain_noDupes[i], format='pdf')
    plt.close()
    np.savetxt("%s overall_biomass.csv" % strain_noDupes[i], zip(binmeans, snow_biomass, floc_biomass), delimiter=",", header="binmeans,%totalbiomass")

