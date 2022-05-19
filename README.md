Requirments:
- Homebrew 2.7.5
- clustalw 2.1
- barrnap 0.9
- bedtools 2.29.2
- tRNAscan-SE 2.0.7
- R 4.0.3
- Python 3.8
- Parallel 20161222

Python3-libraries:
- pandas
- xlsxwriter

R-packages:
- phytools
- xlsx
- rJava
- dplyr
- ggplot2
- gganinimate
- plotly
- vegan
- viridis
- gplots
- heatmaply

De R-packages kunnen geinstalleerd worden via:
install.packages(c('phytools', 'xlsx', 'rJava', 'dplyr', 'ggplot2', 'ggaminimate', 'plotly', 'vegan', 'viridis', 'gplots', 'heatmaply'), dependencies=TRUE, repos = 'http://cran.us.r-project.org')

Voor de pipeline gerunt kan worden moet de input worden aangegeven.
De input moet worden aangegeven in de configfile -> config.yaml.
Inputs die worden gevraagd van de gebruiker zijn:
files -> hier komen de file namen van de genomen die de gebruiker wilt onderzoeken. Dit gaan kom .fna bestanden. 
soort -> hierin moet worden aangegeven of de input Archaea (A), Eukaryoot (E) of Bacterien (B) is. Als dit niet duidelijk is kan de default worden aangehouden; E.
soort_naam -> hier moet de soort input opnieuw worden aangegeven, maar dan voluit geschreven (archaea, bacteria of eukaryote)
Dit moet op de volgende manier worden neergezet:

files:
  - input_file.fna
  - input_file2.fna

GtRNAdbsets:
  - archaea_GtRNAdb
  - bacteria_GtRNAdb
  - eukaryota_GtRNAdb
  
soort:
  - A
  
soort_naam:
  - archaea

LET OP: GtRNAdbsets moet onveranderd blijven!!

De pipeline kan worden aageroepen met de command:
snakemake all

Wanneer alleen de compleetheid moet worden berekend (zonder boom):
snakemake completeness

Wanneer alleen de boom moet worden gemaakt (zonder completeness):
snakemake tree

NMDS-plot: (LET OP: deze wordt niet meegenomen in rule all)
snakemake NMDSwithsubmitted

Rules kunnen ook apart van elkaar gerunt worden. Dat kan met de commando:
snakemake rule_naam

De verschillende rules die worden aangeroepen met all:

- GC_percentage
Deze rule bepaald het GC percentage van de input files van de gebruiker.

- tRNA_scan
Deze rule voert een tRNA scan uit op de input van de gebruiker. Het programma wat hiervoor wordt gebruikt is tRNAscan-SE 2.0.7.

- data_filteren
De input van de gebruiker wordt gefilterd. Alles wat niet uit een A, G, C of T bestaat wordt verwijderd.

- R_op_filtered_data
Fisher exact test wordt uitgevoerd op de gefilterde input van de gebruiker

- fisher_naar_excel
De fisher test resultaten worden weggeschreven naar een excel bestand.

- cosine_simularity
In deze rule wordt de cosine simularity berekend van de gefilterde input van de gebruiker.

- hamming_distance
In deze rule wordt de hamming distance berekend van de gefilterde input van de gebruiker.

- GtRNA_fileren
De tRNA database bestanden van de drie domeinen worden gefilterd.

- NMDS_plot
De NMDS plot wordt gemaakt van alle drie de domeinen.

- histogram_anticodons
Een histogram wordt gemaakt van alle verwachtte anticodons in het domein naar interesse.

- Heatmaps
Er worden heatmaps gemaakt van de anticodons van het domein naar interesse.

- completeness
De completeness wordt geschat aan de hand van de verwachtte anticodons in het volledige domein naar interesse.

- ncbi_genomes
De complete genomen van het domain naar interesse worden opgehaald van NCBI

- genomes_tRNAscan
De complete genomen worden samen met de input van de gebruiker nogmaals door de tRNAscan-SE gehaald.

- boom_genereren
De bestanden voor het genereren van de fylogenetische boom worden aangemaakt. Er wordt met clustalw een multiple sequence alignment gedaan (MSA).

- safe_data
Hier worden gestanden aangemaakt die aangeven in welke genoom welke anticodons aanwezig zijn of juist afwezig.

- tree
De fylogenetische boom wordt gemaakt. In deze boom is te zien waar de input van de gebruiker staan en welke anticodons waar aanwezig of afwezig zijn.


Rules buiten rule all. Deze rules kunnen apart uitgevoerd worden als de gebruiker deze resultaten graag zou willen hebben.
Overige rules:

- NMDSwithsubmitted
Deze rule maakt een NMDS plot aan die de input van de gebruiker gebruikt en de drie domeinen.


OUTPUTS:
Alle bestanden die worden aangemaakt tijdens het aanroepen van de pipeline komen in de directory van de snakefile te staan.