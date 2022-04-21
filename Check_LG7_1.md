# Check the repeats at the end on the 5th contig of chromosome 7 and compare to the ends of the neighboring contigs

## Extract LG7 contig 5, the end of LG7 contig 4 and the start of LG7 contig 6
```
/bin/Python/extractOneSeq.py -i ../GCA_003314205.2_INRA_AMelMel_1.1_genomic.fna -n CM010325.2 -s 4257731 -e 5110297 -o LG7_4
/bin/Python/extractOneSeq.py -i ../GCA_003314205.2_INRA_AMelMel_1.1_genomic.fna -n CM010325.2 -s 5110396 -e 7501794 -o LG7_5
/bin/Python/extractOneSeq.py -i ../GCA_003314205.2_INRA_AMelMel_1.1_genomic.fna -n CM010325.2 -s 7501893 -e 8501893 -o LG7_6
```

## Run Tandem Repeat Finder

```
trf LG7_4 2 7 7 80 10 50 2000 -d -h
trf LG7_5 2 7 7 80 10 50 2000 -d -h
trf LG7_6 2 7 7 80 10 50 2000 -d -h
```



```python
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import re
import csv
import rpy2
from collections import defaultdict
from IPython.display import display

```


```python
dataDict = defaultdict(list)
with open('/Users/avignal/GenotoulSave/Genomes/Abeille/Apis_mellifera_AMelMel1_1/check_LG7/LG7_4.2.7.7.80.10.50.2000.dat') as csvFile:
	data=csv.reader(csvFile, delimiter = " ")
	for row in data:
		if row:
			if re.search('Sequence',row[0]):
				chromosome = row[1]
			elif re.search('\d',row[0]): # and re.match('CM', chromosome) :
				dataDict["GenbankName"].append(chromosome)
				dataDict["RepStart"].append(int(row[0]))
				dataDict["RepEnd"].append(int(row[1]))
				dataDict["Period_Size"].append(int(row[2]))
				dataDict["Copy_Number"].append(float(row[3]))
				dataDict["Percent_Matches"].append(int(row[5]))
				dataDict["Percent_Indels"].append(int(row[6]))
				dataDict["Repseq"].append(row[13])
				dataDict["StartMb"].append(int(row[0])/1000000)
				dataDict["EndMb"].append(int(row[1])/1000000)
				dataDict["MeanPosMb"].append((int(row[0])+int(row[1]))/2000000)
contig4 = pd.DataFrame.from_dict(dataDict, )
contig4['Contig'] = "chr7_contig4"
contig4.loc[:,'fasta'] = ">" + contig4.GenbankName +"_" + contig4.RepStart.map(str) + "_" + contig4.RepEnd.map(str) + "\n" + contig4.Repseq

dataDict = defaultdict(list)
with open('/Users/avignal/GenotoulSave/Genomes/Abeille/Apis_mellifera_AMelMel1_1/check_LG7/LG7_5.2.7.7.80.10.50.2000.dat') as csvFile:
	data=csv.reader(csvFile, delimiter = " ")
	for row in data:
		if row:
			if re.search('Sequence',row[0]):
				chromosome = row[1]
			elif re.search('\d',row[0]): # and re.match('CM', chromosome) :
				dataDict["GenbankName"].append(chromosome)
				dataDict["RepStart"].append(int(row[0]))
				dataDict["RepEnd"].append(int(row[1]))
				dataDict["Period_Size"].append(int(row[2]))
				dataDict["Copy_Number"].append(float(row[3]))
				dataDict["Percent_Matches"].append(int(row[5]))
				dataDict["Percent_Indels"].append(int(row[6]))
				dataDict["Repseq"].append(row[13])
				dataDict["StartMb"].append(int(row[0])/1000000)
				dataDict["EndMb"].append(int(row[1])/1000000)
				dataDict["MeanPosMb"].append((int(row[0])+int(row[1]))/2000000)
contig5 = pd.DataFrame.from_dict(dataDict, )
contig5['Contig'] = "chr7_contig5"
contig5.loc[:,'fasta'] = ">" + contig5.GenbankName +"_" + contig5.RepStart.map(str) + "_" + contig5.RepEnd.map(str) + "\n" + contig5.Repseq

dataDict = defaultdict(list)
with open('/Users/avignal/GenotoulSave/Genomes/Abeille/Apis_mellifera_AMelMel1_1/check_LG7/LG7_6.2.7.7.80.10.50.2000.dat') as csvFile:
	data=csv.reader(csvFile, delimiter = " ")
	for row in data:
		if row:
			if re.search('Sequence',row[0]):
				chromosome = row[1]
			elif re.search('\d',row[0]): # and re.match('CM', chromosome) :
				dataDict["GenbankName"].append(chromosome)
				dataDict["RepStart"].append(int(row[0]))
				dataDict["RepEnd"].append(int(row[1]))
				dataDict["Period_Size"].append(int(row[2]))
				dataDict["Copy_Number"].append(float(row[3]))
				dataDict["Percent_Matches"].append(int(row[5]))
				dataDict["Percent_Indels"].append(int(row[6]))
				dataDict["Repseq"].append(row[13])
				dataDict["StartMb"].append(int(row[0])/1000000)
				dataDict["EndMb"].append(int(row[1])/1000000)
				dataDict["MeanPosMb"].append((int(row[0])+int(row[1]))/2000000)
contig6 = pd.DataFrame.from_dict(dataDict, )
contig6['Contig'] = "chr7_contig6"
contig6.loc[:,'fasta'] = ">" + contig6.GenbankName +"_" + contig6.RepStart.map(str) + "_" + contig6.RepEnd.map(str) + "\n" + contig6.Repseq

```

## Compare end of contig 4 with start of contig 5


```python
contig5[contig5.Period_Size > 50].sort_values(by='RepEnd',ascending=True).head(10)

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GenbankName</th>
      <th>RepStart</th>
      <th>RepEnd</th>
      <th>Period_Size</th>
      <th>Copy_Number</th>
      <th>Percent_Matches</th>
      <th>Percent_Indels</th>
      <th>Repseq</th>
      <th>StartMb</th>
      <th>EndMb</th>
      <th>MeanPosMb</th>
      <th>Contig</th>
      <th>fasta</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>0</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>2</td>
      <td>3563</td>
      <td>258</td>
      <td>13.8</td>
      <td>98</td>
      <td>0</td>
      <td>TTTCGAAATCAGTACGAGCGACGATTTCTTCCGAAAAATTCAATAA...</td>
      <td>0.000002</td>
      <td>0.003563</td>
      <td>0.001783</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_2_3563\nTTTCGAAATC...</td>
    </tr>
    <tr>
      <td>1</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>2</td>
      <td>3563</td>
      <td>515</td>
      <td>6.9</td>
      <td>98</td>
      <td>0</td>
      <td>TTTCGAAATCAGTACGATCGACGATTTCTTCCGAAAAATTCAATAA...</td>
      <td>0.000002</td>
      <td>0.003563</td>
      <td>0.001783</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_2_3563\nTTTCGAAATC...</td>
    </tr>
    <tr>
      <td>2</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>3346</td>
      <td>3784</td>
      <td>221</td>
      <td>2.0</td>
      <td>98</td>
      <td>0</td>
      <td>CTTTCGAAATCAGTACGAGCGACGATTTCTTCCGAAAAATTCAATA...</td>
      <td>0.003346</td>
      <td>0.003784</td>
      <td>0.003565</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_3346_3784\nCTTTCGA...</td>
    </tr>
    <tr>
      <td>3</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>3567</td>
      <td>17103</td>
      <td>258</td>
      <td>52.5</td>
      <td>98</td>
      <td>0</td>
      <td>CTTTCGAAATCAGTACGAGCGACGATTTCTTCCGAAAAATTCAATA...</td>
      <td>0.003567</td>
      <td>0.017103</td>
      <td>0.010335</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_3567_17103\nCTTTCG...</td>
    </tr>
    <tr>
      <td>158</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>84449</td>
      <td>84673</td>
      <td>109</td>
      <td>2.1</td>
      <td>100</td>
      <td>0</td>
      <td>AAAATAACATGTAAACTTATAAATACAAACATAAATTGAATGCAAT...</td>
      <td>0.084449</td>
      <td>0.084673</td>
      <td>0.084561</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_84449_84673\nAAAAT...</td>
    </tr>
    <tr>
      <td>216</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>106139</td>
      <td>106289</td>
      <td>56</td>
      <td>2.8</td>
      <td>86</td>
      <td>10</td>
      <td>ATTATAATTATATAATTTATATAATTATATAATTATTATTATTATA...</td>
      <td>0.106139</td>
      <td>0.106289</td>
      <td>0.106214</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_106139_106289\nATT...</td>
    </tr>
    <tr>
      <td>376</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>178251</td>
      <td>178437</td>
      <td>88</td>
      <td>2.1</td>
      <td>94</td>
      <td>5</td>
      <td>AAAATATAAGCTAGCATTATCTTTAAAATCAAATAAAGCTAATTTT...</td>
      <td>0.178251</td>
      <td>0.178437</td>
      <td>0.178344</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_178251_178437\nAAA...</td>
    </tr>
    <tr>
      <td>957</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>456571</td>
      <td>456715</td>
      <td>64</td>
      <td>2.3</td>
      <td>90</td>
      <td>9</td>
      <td>ATATTATTTAATATATATATAATTTAATATATATATATATATTATT...</td>
      <td>0.456571</td>
      <td>0.456715</td>
      <td>0.456643</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_456571_456715\nATA...</td>
    </tr>
    <tr>
      <td>1154</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>541517</td>
      <td>541685</td>
      <td>83</td>
      <td>2.0</td>
      <td>100</td>
      <td>0</td>
      <td>AATTTGATGAAATAATACATATTTCAATATTTCTTTAAAATAAAGT...</td>
      <td>0.541517</td>
      <td>0.541685</td>
      <td>0.541601</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_541517_541685\nAAT...</td>
    </tr>
    <tr>
      <td>1307</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>606514</td>
      <td>606795</td>
      <td>141</td>
      <td>2.0</td>
      <td>100</td>
      <td>0</td>
      <td>AAGTTCTCTAAATACTTCTTTCAATTCTATTCTTCTTTCATTAAAT...</td>
      <td>0.606514</td>
      <td>0.606795</td>
      <td>0.606654</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_606514_606795\nAAG...</td>
    </tr>
  </tbody>
</table>
</div>




```python
contig4[contig4.Period_Size > 50].sort_values(by='RepEnd',ascending=False).head(10)

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GenbankName</th>
      <th>RepStart</th>
      <th>RepEnd</th>
      <th>Period_Size</th>
      <th>Copy_Number</th>
      <th>Percent_Matches</th>
      <th>Percent_Indels</th>
      <th>Repseq</th>
      <th>StartMb</th>
      <th>EndMb</th>
      <th>MeanPosMb</th>
      <th>Contig</th>
      <th>fasta</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1738</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>838236</td>
      <td>852566</td>
      <td>773</td>
      <td>18.5</td>
      <td>98</td>
      <td>0</td>
      <td>GTTTGGCAGTATCACAATTTCGTATGTATAAGACGCTTGATTTCAT...</td>
      <td>0.838236</td>
      <td>0.852566</td>
      <td>0.845401</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_838236_852566\nGTT...</td>
    </tr>
    <tr>
      <td>1737</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>838236</td>
      <td>852566</td>
      <td>258</td>
      <td>55.6</td>
      <td>98</td>
      <td>0</td>
      <td>GTTTGGCAGTATCACAATTTCGTATGTATAAGACGCTTGATTTCAT...</td>
      <td>0.838236</td>
      <td>0.852566</td>
      <td>0.845401</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_838236_852566\nGTT...</td>
    </tr>
    <tr>
      <td>1727</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>833762</td>
      <td>833915</td>
      <td>53</td>
      <td>2.7</td>
      <td>85</td>
      <td>12</td>
      <td>TATTATTATTATATTATATTATATTTATATATTATATATATATATA...</td>
      <td>0.833762</td>
      <td>0.833915</td>
      <td>0.833839</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_833762_833915\nTAT...</td>
    </tr>
    <tr>
      <td>1705</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>822694</td>
      <td>822798</td>
      <td>54</td>
      <td>1.9</td>
      <td>92</td>
      <td>0</td>
      <td>TGACGATCTCGACTAATATGTCTATCACGACTAATATGTCGATCTC...</td>
      <td>0.822694</td>
      <td>0.822798</td>
      <td>0.822746</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_822694_822798\nTGA...</td>
    </tr>
    <tr>
      <td>1451</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>685323</td>
      <td>685510</td>
      <td>82</td>
      <td>2.3</td>
      <td>83</td>
      <td>11</td>
      <td>ATATATATTATATATATTATTATAATATATATATTTATTTATATAT...</td>
      <td>0.685323</td>
      <td>0.685510</td>
      <td>0.685416</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_685323_685510\nATA...</td>
    </tr>
    <tr>
      <td>1286</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>619846</td>
      <td>620142</td>
      <td>152</td>
      <td>2.0</td>
      <td>92</td>
      <td>7</td>
      <td>ATATTATTCTTCTATATTAAAAAAAATATTTTTTTTATCATGATCA...</td>
      <td>0.619846</td>
      <td>0.620142</td>
      <td>0.619994</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_619846_620142\nATA...</td>
    </tr>
    <tr>
      <td>1090</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>540294</td>
      <td>540588</td>
      <td>159</td>
      <td>1.8</td>
      <td>88</td>
      <td>1</td>
      <td>ACTAAACACATCTAATTAATAAAAAAGATTTAGGTTAAGATGCTAT...</td>
      <td>0.540294</td>
      <td>0.540588</td>
      <td>0.540441</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_540294_540588\nACT...</td>
    </tr>
    <tr>
      <td>1080</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>535590</td>
      <td>535704</td>
      <td>57</td>
      <td>2.0</td>
      <td>100</td>
      <td>0</td>
      <td>TTTTTGTATTTTTTTTAATATAGTTAAATATATGAAATCATAAAAT...</td>
      <td>0.535590</td>
      <td>0.535704</td>
      <td>0.535647</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_535590_535704\nTTT...</td>
    </tr>
    <tr>
      <td>978</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>481489</td>
      <td>481600</td>
      <td>56</td>
      <td>2.0</td>
      <td>98</td>
      <td>1</td>
      <td>AATTATTTCTGAGTATAATATTATAATTAATATTATATATTATAAT...</td>
      <td>0.481489</td>
      <td>0.481600</td>
      <td>0.481544</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_481489_481600\nAAT...</td>
    </tr>
    <tr>
      <td>971</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>479831</td>
      <td>479957</td>
      <td>55</td>
      <td>2.3</td>
      <td>91</td>
      <td>5</td>
      <td>TAAATTAAATTTTTATTTTTTTTTTATATTAAAATCATAATATTTA...</td>
      <td>0.479831</td>
      <td>0.479957</td>
      <td>0.479894</td>
      <td>chr7_contig4</td>
      <td>&gt;CM010325.2_4257731_5110297_479831_479957\nTAA...</td>
    </tr>
  </tbody>
</table>
</div>



* The repeat of period size 773 at the end of contig 4 has the same coordinates as the repeat of period size 258 and 258 * 3 = 774.
* Pairwise alignment of the 258 repeat on the 773 repeat gives 3 possible contigous alignments with 100%, 100% and 99.61% (1 mismatch) identity.
* They are thus equivalent.
## Extract sequence of the repeat found at the ends of both contigs, of period size 258


```python
fasta1 = contig5.loc[(contig5.RepStart == 2) & (contig5.RepEnd == 3563) & (contig5.Period_Size == 258)]
fasta2 = contig5.loc[(contig5.RepStart == 3567) & (contig5.RepEnd == 17103) & (contig5.Period_Size == 258)]
fasta3 = contig4.loc[(contig4.RepStart == 838236) & (contig4.RepEnd == 852566) & (contig4.Period_Size == 258)]

```

## Compare end of contig 5 with start of contig 6


```python
contig5[contig5.Period_Size > 50].sort_values(by='RepEnd',ascending=False).head(10)

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GenbankName</th>
      <th>RepStart</th>
      <th>RepEnd</th>
      <th>Period_Size</th>
      <th>Copy_Number</th>
      <th>Percent_Matches</th>
      <th>Percent_Indels</th>
      <th>Repseq</th>
      <th>StartMb</th>
      <th>EndMb</th>
      <th>MeanPosMb</th>
      <th>Contig</th>
      <th>fasta</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>4225</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>2375287</td>
      <td>2391398</td>
      <td>1296</td>
      <td>12.5</td>
      <td>98</td>
      <td>1</td>
      <td>AATTCATCATTTCATTCTATTTTTCAATCTCGAAATACATTATGGC...</td>
      <td>2.375287</td>
      <td>2.391398</td>
      <td>2.383342</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_2375287_2391398\nA...</td>
    </tr>
    <tr>
      <td>4153</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>2259517</td>
      <td>2259842</td>
      <td>164</td>
      <td>2.0</td>
      <td>95</td>
      <td>0</td>
      <td>CTCGACGCATCCTGCCTTTTGTTCGCCGATAACCAAGAATTAATAA...</td>
      <td>2.259517</td>
      <td>2.259842</td>
      <td>2.259679</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_2259517_2259842\nC...</td>
    </tr>
    <tr>
      <td>4077</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>2175140</td>
      <td>2175254</td>
      <td>62</td>
      <td>1.9</td>
      <td>87</td>
      <td>8</td>
      <td>TATATTATATTAATATAATATATTATTATTATATATTATTATATAT...</td>
      <td>2.175140</td>
      <td>2.175254</td>
      <td>2.175197</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_2175140_2175254\nT...</td>
    </tr>
    <tr>
      <td>3959</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>2088582</td>
      <td>2088745</td>
      <td>79</td>
      <td>2.1</td>
      <td>100</td>
      <td>0</td>
      <td>AAGAAAATACTGAATTAAATTAAGTAATTATGAATATGAAAAAATG...</td>
      <td>2.088582</td>
      <td>2.088745</td>
      <td>2.088664</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_2088582_2088745\nA...</td>
    </tr>
    <tr>
      <td>3735</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>1908341</td>
      <td>1908449</td>
      <td>52</td>
      <td>2.1</td>
      <td>84</td>
      <td>6</td>
      <td>ATTATTTATTATTAATATTAATAATATTATTAATATTAATAATATT...</td>
      <td>1.908341</td>
      <td>1.908449</td>
      <td>1.908395</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_1908341_1908449\nA...</td>
    </tr>
    <tr>
      <td>3606</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>1810249</td>
      <td>1810437</td>
      <td>99</td>
      <td>1.9</td>
      <td>86</td>
      <td>6</td>
      <td>ATTTATAATATATTATATATATATAAATATATTATATATATATATA...</td>
      <td>1.810249</td>
      <td>1.810437</td>
      <td>1.810343</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_1810249_1810437\nA...</td>
    </tr>
    <tr>
      <td>3533</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>1770992</td>
      <td>1771108</td>
      <td>52</td>
      <td>2.3</td>
      <td>95</td>
      <td>4</td>
      <td>ATATAATATAATATACATATATTAAATTGTTAATTTAATTTTATAT...</td>
      <td>1.770992</td>
      <td>1.771108</td>
      <td>1.771050</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_1770992_1771108\nA...</td>
    </tr>
    <tr>
      <td>3272</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>1660876</td>
      <td>1661117</td>
      <td>120</td>
      <td>2.0</td>
      <td>99</td>
      <td>0</td>
      <td>AATTGATCTTGATTTTAAAATTGAATTAATCAAATAAATATTTAAT...</td>
      <td>1.660876</td>
      <td>1.661117</td>
      <td>1.660996</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_1660876_1661117\nA...</td>
    </tr>
    <tr>
      <td>3249</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>1640486</td>
      <td>1640895</td>
      <td>111</td>
      <td>3.7</td>
      <td>86</td>
      <td>7</td>
      <td>ATATATATATATATATATACTATTATATATACATACATATATATAT...</td>
      <td>1.640486</td>
      <td>1.640895</td>
      <td>1.640691</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_1640486_1640895\nA...</td>
    </tr>
    <tr>
      <td>3224</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>1632077</td>
      <td>1632290</td>
      <td>112</td>
      <td>1.9</td>
      <td>100</td>
      <td>0</td>
      <td>ATATATATATATATATATATATATATATATATATATATATATATAT...</td>
      <td>1.632077</td>
      <td>1.632290</td>
      <td>1.632184</td>
      <td>chr7_contig5</td>
      <td>&gt;CM010325.2_5110396_7501794_1632077_1632290\nA...</td>
    </tr>
  </tbody>
</table>
</div>




```python
contig6[contig6.Period_Size > 50].sort_values(by=['RepStart']).head(10)

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GenbankName</th>
      <th>RepStart</th>
      <th>RepEnd</th>
      <th>Period_Size</th>
      <th>Copy_Number</th>
      <th>Percent_Matches</th>
      <th>Percent_Indels</th>
      <th>Repseq</th>
      <th>StartMb</th>
      <th>EndMb</th>
      <th>MeanPosMb</th>
      <th>Contig</th>
      <th>fasta</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>0</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>4</td>
      <td>4649</td>
      <td>1296</td>
      <td>3.6</td>
      <td>99</td>
      <td>0</td>
      <td>ATAACAAAATATAAATACTGCTACTTATTATCATATATTATTCAAG...</td>
      <td>0.000004</td>
      <td>0.004649</td>
      <td>0.002327</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_4_4649\nATAACAAAAT...</td>
    </tr>
    <tr>
      <td>40</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>73111</td>
      <td>74373</td>
      <td>306</td>
      <td>4.2</td>
      <td>93</td>
      <td>3</td>
      <td>ATTAATGGACTATTAAAAAAAGGAACATTAAAACAGTTGAAAAAAA...</td>
      <td>0.073111</td>
      <td>0.074373</td>
      <td>0.073742</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_73111_74373\nATTAA...</td>
    </tr>
    <tr>
      <td>41</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>73173</td>
      <td>74373</td>
      <td>613</td>
      <td>2.0</td>
      <td>91</td>
      <td>4</td>
      <td>AAAAACAACTATCAGTCAAAATCTATTTGTCTATGGAATAAATATA...</td>
      <td>0.073173</td>
      <td>0.074373</td>
      <td>0.073773</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_73173_74373\nAAAAA...</td>
    </tr>
    <tr>
      <td>83</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>125600</td>
      <td>126016</td>
      <td>224</td>
      <td>1.9</td>
      <td>91</td>
      <td>5</td>
      <td>TAAATATCGATATAACGATTATACTTTAAGAGTTTGAAGAGAAAAA...</td>
      <td>0.125600</td>
      <td>0.126016</td>
      <td>0.125808</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_125600_126016\nTAA...</td>
    </tr>
    <tr>
      <td>99</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>158731</td>
      <td>158853</td>
      <td>60</td>
      <td>2.0</td>
      <td>98</td>
      <td>1</td>
      <td>AGAAATAATTTTACAGAATATGTATATATTTAACATGTTTGAGAAA...</td>
      <td>0.158731</td>
      <td>0.158853</td>
      <td>0.158792</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_158731_158853\nAGA...</td>
    </tr>
    <tr>
      <td>102</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>165027</td>
      <td>165830</td>
      <td>310</td>
      <td>2.6</td>
      <td>88</td>
      <td>4</td>
      <td>TATGATCATCACATTTCTCTTTTCAAAAATTAATAAAATTAGTATT...</td>
      <td>0.165027</td>
      <td>0.165830</td>
      <td>0.165429</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_165027_165830\nTAT...</td>
    </tr>
    <tr>
      <td>109</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>170906</td>
      <td>171216</td>
      <td>112</td>
      <td>2.9</td>
      <td>84</td>
      <td>6</td>
      <td>AAAAAATTTATTTACTATTACAAAATTCCAGATTATACTCTGATTC...</td>
      <td>0.170906</td>
      <td>0.171216</td>
      <td>0.171061</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_170906_171216\nAAA...</td>
    </tr>
    <tr>
      <td>122</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>185529</td>
      <td>186733</td>
      <td>299</td>
      <td>4.0</td>
      <td>90</td>
      <td>4</td>
      <td>AAATTTGTTAGGAAGGGTCATAAAATGGAAAGATTTTCGTTCTTTG...</td>
      <td>0.185529</td>
      <td>0.186733</td>
      <td>0.186131</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_185529_186733\nAAA...</td>
    </tr>
    <tr>
      <td>123</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>185529</td>
      <td>186733</td>
      <td>607</td>
      <td>2.0</td>
      <td>89</td>
      <td>6</td>
      <td>AAATTTGTTAGGAAGGGTCATAAAATGGAAAGATTTTCGTTCTTTG...</td>
      <td>0.185529</td>
      <td>0.186733</td>
      <td>0.186131</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_185529_186733\nAAA...</td>
    </tr>
    <tr>
      <td>124</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>190475</td>
      <td>190712</td>
      <td>130</td>
      <td>1.8</td>
      <td>93</td>
      <td>1</td>
      <td>AATTGATATCCAATTTACAATATTATAATTATTTCTATTTTACATG...</td>
      <td>0.190475</td>
      <td>0.190712</td>
      <td>0.190593</td>
      <td>chr7_contig6</td>
      <td>&gt;CM010325.2_7501893_8501893_190475_190712\nAAT...</td>
    </tr>
  </tbody>
</table>
</div>



## Extract sequence of the repeat found at the ends of both contigs, of period size 1296


```python
fasta4 = contig5.loc[(contig5.RepStart == 2375287) & (contig5.RepEnd == 2391398) & (contig5.Period_Size == 1296)]
fasta5 = contig6.loc[(contig6.RepStart == 4) & (contig6.RepEnd == 4649) & (contig6.Period_Size == 1296)]

```

## Export all 5 sequences in fasta format for mafft alignments



```python
all_sequences = pd.concat([fasta1,fasta2,fasta3,fasta4,fasta5], axis = 0)
#all_sequences.fasta.to_csv('/Users/avignal/GenotoulSave/Genomes/Abeille/Apis_mellifera_AMelMel1_1/check_LG7/junctions_1_2.fa', index=False, quoting = csv.QUOTE_NONE, escapechar = ' ', header=False)

```

## Table for supplementary


```python
export = all_sequences[['GenbankName','Contig','RepStart','RepEnd','Period_Size','Copy_Number','Percent_Matches','Percent_Indels','Repseq']]
export = export.sort_values(by=['Contig'])
#export.to_csv('/Users/avignal/Documents/Articles/2019_PacificBee/table_ends_LG7.xls', index=False, sep="\t")
export

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GenbankName</th>
      <th>Contig</th>
      <th>RepStart</th>
      <th>RepEnd</th>
      <th>Period_Size</th>
      <th>Copy_Number</th>
      <th>Percent_Matches</th>
      <th>Percent_Indels</th>
      <th>Repseq</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1737</td>
      <td>CM010325.2_4257731_5110297</td>
      <td>chr7_contig4</td>
      <td>838236</td>
      <td>852566</td>
      <td>258</td>
      <td>55.6</td>
      <td>98</td>
      <td>0</td>
      <td>GTTTGGCAGTATCACAATTTCGTATGTATAAGACGCTTGATTTCAT...</td>
    </tr>
    <tr>
      <td>0</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>chr7_contig5</td>
      <td>2</td>
      <td>3563</td>
      <td>258</td>
      <td>13.8</td>
      <td>98</td>
      <td>0</td>
      <td>TTTCGAAATCAGTACGAGCGACGATTTCTTCCGAAAAATTCAATAA...</td>
    </tr>
    <tr>
      <td>3</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>chr7_contig5</td>
      <td>3567</td>
      <td>17103</td>
      <td>258</td>
      <td>52.5</td>
      <td>98</td>
      <td>0</td>
      <td>CTTTCGAAATCAGTACGAGCGACGATTTCTTCCGAAAAATTCAATA...</td>
    </tr>
    <tr>
      <td>4225</td>
      <td>CM010325.2_5110396_7501794</td>
      <td>chr7_contig5</td>
      <td>2375287</td>
      <td>2391398</td>
      <td>1296</td>
      <td>12.5</td>
      <td>98</td>
      <td>1</td>
      <td>AATTCATCATTTCATTCTATTTTTCAATCTCGAAATACATTATGGC...</td>
    </tr>
    <tr>
      <td>0</td>
      <td>CM010325.2_7501893_8501893</td>
      <td>chr7_contig6</td>
      <td>4</td>
      <td>4649</td>
      <td>1296</td>
      <td>3.6</td>
      <td>99</td>
      <td>0</td>
      <td>ATAACAAAATATAAATACTGCTACTTATTATCATATATTATTCAAG...</td>
    </tr>
  </tbody>
</table>
</div>




```python

```
