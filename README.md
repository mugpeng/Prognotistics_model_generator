# Updates

## 0327

fix a function bug of wilcoxDEGs.



new:

```
wilcoxDEGs <- function(exp_df, clinical1){
  # exp_df in log2 format
  exp_df <- 2**(exp_df) - 1
  re <- lapply(1:nrow(exp_df), function(x){
    id1 <- clinical1[clinical1$group %in% levels(clinical1$group)[1],]$patient
    id2 <- clinical1[clinical1$group %in% levels(clinical1$group)[2],]$patient
    group1 <- exp_df[,colnames(exp_df) %in% id1]
    group2 <- exp_df[,colnames(exp_df) %in% id2]
    p <- wilcox.test(as.numeric(group2[x,]), as.numeric(group1[x,]))$p.value
    logFC <- log2(mean(as.numeric(group2[x,])) / mean(as.numeric(group1[x,])))
    c(p, logFC)
  })
  re <- as.data.frame(do.call(rbind, re))
  rownames(re) <- rownames(exp_df)
  colnames(re) <- c("P.value", "LogFC")
  re$FDR <- p.adjust(re$P.value, method = "fdr") 
  return(re)
}
```



old:

```
wilcoxDEGs <- function(exp_df, clinical1){
  re <- lapply(1:nrow(exp_df), function(x){
    id1 <- clinical1[clinical1$group %in% levels(clinical1$group)[1],]$patient
    id2 <- clinical1[clinical1$group %in% levels(clinical1$group)[2],]$patient
    group1 <- exp_df[,colnames(exp_df) %in% id1]
    group2 <- exp_df[,colnames(exp_df) %in% id2]
    p <- wilcox.test(as.numeric(group2[x,]), as.numeric(group1[x,]))$p.value
    logFC <- mean(as.numeric(group2[x,])) - mean(as.numeric(group1[x,]))
    c(p, logFC)
  })
  re <- as.data.frame(do.call(rbind, re))
  rownames(re) <- rownames(exp_df)
  colnames(re) <- c("P.value", "LogFC")
  re$FDR <- p.adjust(re$P.value, method = "fdr") 
  return(re)
}
```



The change may influence some models outputs, but my wrong wilcox function still capture the large exp diff genes. 



# Some mumbles

It's just a simple attempt that three very simple ML methods with different parameters, input selection(tumor types, feature selection) can generate 118,956,840 models. And there are numerous prognostics models created for publishing without any validations or clinical applications.

[There is no such thing as a validated prediction model | BMC Medicine | Full Text](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-023-02779-w#Sec2)

What if I generate many many many models? 



# Design

Here, I create TPMG(tumor prognostics model generator).

The pan RNA-seq data including 38 tumor types were collected from xena (https://xenabrowser.net/), the Fragments Per Kilobase of transcript per Million mapped reads upper quartile (FPKM-UQ) with log-2 transformed was used for normalization. The RNA-seq data was processed to remain mRNA expression by filtering genes with "protein_coding" in gencode annotation gtf file(v41).

Survival data corresponding with the RNA-seq data was obtained from xena (https://xenabrowser.net/) and TARGET (https://target-data.nci.nih.gov/Public/), we kept patients with overall survival(OS) and overall survival time(OS.time) data.

Gene pathways were downloaded from  The Molecular Signatures Database (MSigDB) of Human Collections (https://www.gsea-msigdb.org/gsea/msigdb). For selecting meaningful pathways, we retained pathways from HALLMARK(50), KEGG(186), GO_BP(7763) and GO_CC(1035), GO_MF(1763), which including 10,797 pathways.

For developing a satisfied prognostic model, I performed a 3-fold CV with 30 repeats of three machine learning methods including, lasso regression model(R package "glmnet"), cox regression(R package "survival") and random forest(R package "randomForest"). I also applied these models with different parameters. Specifically, in lasso regression model, the parameter is alpha = 0,0.1,0.2,...,0.8,0.9,1; in random forest. the parameter is mtry = 10, 15, 20, 30, ntree = 20, 50, 100, 200; in cox regression model, the parameter is step wise step = "yes" or "no". So theoretically, It can produce 870(30 × 11 + 30 × 4 × 4 + 30 × 2) models for each pathway, 9,393,390(870 × 10,797) models for each tumor, 356,948,820(9,393,390 × 38) models altogether.



# Some descriptions

- Theoretically, the TPMG can produce 347,555,430 models for each tumor-pathway pairs. But some pathways might not have sufficient genes to generate models or less sufficient thus been discarded.
- TGCT tumor is missing because too less survival records to generate models.
- You can change the scripts for your tumor data or interested pathways(gene list).
- The input and temporary files, model generator results(parameter) are uploaded on figshare.



# Code orders

- Preparation files:
  - Function.R
  - bash_script, for create files for multi-threads



- Data preprocess
  - R01
- Model generators 
  - Main.R, iterate all pathways in different tumors；
  - bash_script/make_parallel_script.sh，；
  - bash_script/split_all_script.sh，
    - -l for assigning how many pathways to run in a script，(total pathway counts/i) is the parallel number；
  - script_run.sh，run these script when everything is ready；
- Validate model from outside
  - only for lasso validation right now，R02
- Subsequent analysis including survival, immune infiltration analysis.
  - R03、R04
- Statics analysis 
  - R05、R06、R07
  - Are there any possibilities that the best performance pathways or the intersected genes in these pathways can reveal some new findings like cancerous genes?



# Results

## Basic statistics

I finally obtained 118,956,840 models, and these models were splited by tumor types and uploaded at [Model_generators_results](https://doi.org/10.6084/m9.figshare.22247455.v5). 

We calculate the median AUC of each method, lasso model(AUC, 0.682) has the highest AUC compared to cox(AUC, 0.650) and random forest(AUC, 0.610).

![](http://cos01.mugpeng.top/img/f1.png)

Split by tumor types:

![](http://cos01.mugpeng.top/img/f2.png)



Adrenocortical Cancer(ACC), Kidney Chromophobe(KICH) and Ocular melanomas(UVM) had highest AUC, while those from Lung Squamous Cell Carcinoma(LUSC), Stomach Cancer(STAD), Breast Cancer(BRCA) and Glioblastoma(GBM) had lowest AUC:

![](http://cos01.mugpeng.top/img/f3.png)



From my personal experience, AUC > 0.7 could be a satisfied threshold, in this case, I remained models with  AUC larger than 0.7 both in testing and training. After filtering 23,849,573 models retained, and there was no models for LUSC. 

![](http://cos01.mugpeng.top/img/f4.png)



## Pathways

Next, I used the remained 23,849,573 models and extracted the top20 median AUC pathways in each tumor, the result indicated that the highest frequent pathways including GOBP_ION_TRANSPORT, GOMF_TRANSPORTER_ACTIVITY, GOCC_MEMBRANE_PROTEIN_COMPLEX...:

```

                    GOBP_SYNAPTIC_SIGNALING 
                                          4 
                GOCC_PLASMA_MEMBRANE_REGION 
                                          4 
                               GOCC_SYNAPSE 
                                          4 
            GOMF_SIGNALING_RECEPTOR_BINDING 
                                          4 
GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND 
                                          5 
               GOBP_TRANSMEMBRANE_TRANSPORT 
                                          5 
              GOCC_MEMBRANE_PROTEIN_COMPLEX 
                                          5 
                  GOMF_TRANSPORTER_ACTIVITY 
                                          5 
           GOBP_ION_TRANSMEMBRANE_TRANSPORT 
                                          7 
                         GOBP_ION_TRANSPORT 
                                          7 
```



While among the remained 23,849,573 models, there are a huge amount of pathways:

![](http://cos01.mugpeng.top/img/f5.png)

And I selected the most frequency pathway, it appeared that there were 210 pathways at least in 30 tumors, some pathways were:

```

                               GOBP_CELLULAR_RESPONSE_TO_ENDOGENOUS_STIMULUS 
                                                                          33 
                                                  GOBP_GENERATION_OF_NEURONS 
                                                                          33 
                                                             GOBP_LOCOMOTION 
                                                                          33 
                           GOBP_NEGATIVE_REGULATION_OF_DEVELOPMENTAL_PROCESS 
                                                                          33 
                GOBP_NEGATIVE_REGULATION_OF_MULTICELLULAR_ORGANISMAL_PROCESS 
                                                                          33 
GOBP_NEGATIVE_REGULATION_OF_NUCLEOBASE_CONTAINING_COMPOUND_METABOLIC_PROCESS 
                                                                          33 
                           GOBP_NEGATIVE_REGULATION_OF_RNA_METABOLIC_PROCESS 
                                                                          33 
                                                           GOBP_NEUROGENESIS 
                                                                          33 
                           GOBP_POSITIVE_REGULATION_OF_DEVELOPMENTAL_PROCESS 
                                                                          33 
              GOBP_POSITIVE_REGULATION_OF_MACROMOLECULE_BIOSYNTHETIC_PROCESS 
                                                                          33 
                              GOBP_POSITIVE_REGULATION_OF_MOLECULAR_FUNCTION 
                                                                          33 
                GOBP_POSITIVE_REGULATION_OF_MULTICELLULAR_ORGANISMAL_PROCESS 
                                                                          33 
              GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_BY_RNA_POLYMERASE_II 
                                                                          33 
                                     GOBP_REGULATION_OF_CELL_DIFFERENTIATION 
                                                                          33 
                                                      GOBP_RESPONSE_TO_LIPID 
                                                                          33 
                                 GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND 
                                                                          33 
                                                      GOCC_NEURON_PROJECTION 
                                                                          33 
          GOCC_NUCLEAR_OUTER_MEMBRANE_ENDOPLASMIC_RETICULUM_MEMBRANE_NETWORK 
                                                                          33 
                                                          GOMF_LIPID_BINDING 
                                                                          33 
                                          GOMF_SEQUENCE_SPECIFIC_DNA_BINDING 
                                                                          33 
                                             GOMF_SIGNALING_RECEPTOR_BINDING 
                                                                          33 
                                  GOMF_SIGNALING_RECEPTOR_REGULATOR_ACTIVITY 
                                                                          33 
                                                   GOMF_TRANSPORTER_ACTIVITY 
                                                                          33 
                                                          GOBP_CELL_ADHESION 
                                                                          34 
                           GOBP_POSITIVE_REGULATION_OF_RNA_METABOLIC_PROCESS 
                                                                          34 
                                          GOBP_RESPONSE_TO_NITROGEN_COMPOUND 
                                                                          34 
                                                              GOCC_CHROMATIN 
                                                                          34 
```



Interestingly, GOMF_TRANSPORTER_ACTIVITY was present in both "top20" and "all" selection strategy pathway. GOMF_TRANSPORTER_ACTIVITY was consisted of 1231 genes.



Then we selected the genes from GOMF_TRANSPORTER_ACTIVITY and used to construct models of BRCA, CHOL, HNSC, READ, and UCEC. The top frequency genes including:

```
PTEN, WNT5A, SIX1, ADCY8, ADRA2A, IFNG, CLU, TLR9, SOX11, SLC8A1
```





## Genes

the genes of top model of each tumor were listed here:

```
STAD,"ACADL, AK4, ALAS2, ARMCX1, ATP2A1, BDH2, BNIP3, BRI3BP, CAVIN1, CNR1, CRYAB, CYP1B1, ERBB4, FAM110B, FEN1, FEZ1, FKBP10, GABARAPL1, GATD3, GHR, GLYATL2, MAPK10, MICU3, MSRB3, NDRG4, NOX4, NPTX1, NR3C1, PABPC5, PDE2A, PDK4, PLN, RGS2, RPS6KA6, SLC25A21, SLIT3, SNCA, SPARC, TWNK"
BRCA,"ADRB1, ANXA8, APOA5, BCL2L10, CCL1, CST1, CSTA, CXCL13, FGD3, FOXJ1, GPLD1, GSN, KLRC4_KLRK1, KLRK1, LAMP3, LTF, NEDD9, NGFR, NTF3, NTRK2, NTRK3, ODAM, OPRPN, PAPLN, PDE3A, PLIN5, PMAIP1, PPP1R15A, PPP1R17, PPP1R35, PPP4R4, PRKCZ, PSENEN, PSME2, PTGS2, SERPINA12, SERPINA9, SERPINB5, SFRP1, SIAH2, SMR3A, SPINK2, SPINT1, SPP2, TBC1D4, TFPI2, TP63, VAV3, WNT11"
LUAD,"ATP1B2, CTNND2, DRP2, EFNB2, GAP43, GRIA1, GRIK2, GRIP1, GSG1L, NTSR1, PRKCZ, PTPRO, RAPGEF4, RS1, SLITRK1, SORCS2, SRCIN1, TMEM108, ZDHHC15"
LIHC,"ANG, BEND3, CAMK4, CAPG, CDC14B, CDC6, CDCA8, CENPH, CENPW, CLEC3B, CPS1, DTL, EME1, FANCD2, HJURP, HOXD9, MCM10, MCM3, NEK2, NLRP5, NSD2, PANK1, PKMYT1, RCL1, RGS2, SIX1, TCOF1, UBE2T"
OV,"CREB3, ELK3, FOSB, FOXJ2, GCM1, GLIS2, HNF4A, HOXB2, IKZF3, IRF4, JUNB, KLF13, ONECUT1, PHOX2B, RUNX1, TCF15, TFAP2B, TFCP2L1, ZNF485, ZNF574, ZNF76, ZNF780B"
ESCA,"AGRN, AKAP6, BTG2, FFAR2, FMO1, GPER1, GRAMD1C, HCN1, MAOB, PLA2G1B, SLIT3, STC2"
BLCA,"ADAMTS16, ARHGEF4, CAV1, CC2D2A, CDKL1, CFAP53, COBL, CORO1C, CSPG4, DNM3, DOCK11, DRC1, DYNC2H1, DZIP1, EMP1, FGD5, FSCN1, GAP43, HRG, MSTN, RAB23, RAC3, SH3BP1, SLC9A3R1, VCL, WDR90, WWTR1"
HNSC,"ABCB1, ANO1, BEST2, FCRL3, GRIA3, GRIA4, HRG, KCNAB3, NALCN, OSR1, OTOP3, PER2, SCNN1D, SLC20A1, SLC5A9, SLC6A7, SPINK1, STC1, SYT1, TMEM150C, TRPV1, XPR1"
COAD,"ALOX12B, ATP6V1B1, BMP5, CBLIF, CPT2, CYP26A1, ENOPH1, ERFE, FGF19, FGF23, FTCD, G6PC2, GAL3ST3, GSTM1, GSTM2, HPD, HSPA1A, KCNJ11, LEP, LHB, LRP2, OSBPL1A, PHGDH, PLP1, PPARGC1A, PTH1R, RCVRN, RIMKLB, RORC, SLC4A4, SNAI1, SULT1C3, ZNF692"
SKCM,"ADCYAP1R1, ATG5, ATP1A2, ATP1B1, CALM3, CCL2, CHP2, CLIC2, CORO1A, CRHBP, CTSS, CXCL10, CXCL11, CXCL9, DRD3, GLRX, GRM6, HAMP, JSRP1, KCNMB1, KCNN4, MS4A1, PIK3CG, PPIF, PTK2B, SLC31A2, SLC8A1, STAC3, STIM2, TMC2, TMSB4X, UTRN, VDAC1, WNK2, XCL1"
GBM,"ACTR3C, ALKBH4, DBNL, ENAH, EPS8L2, FXYD5, MARCKS, MEFV, MLPH, MYBPC2, MYO15A, NOD2, NOS3, PARVA, PARVB, PDLIM1, PDLIM4, RDX, SSH3, TLNRD1, TNNC2, TNS4, TULP1, VILL"
UCEC,"CDK5R2, CTF1, DGKG, KALRN, NRG3, OLIG3, ONECUT2, OPCML, OSTN, PTEN, RASAL1, SH3GL2, SIX1, SLITRK2, SOX11, STMN4, STXBP1"
PRAD,"ADRA2A, ATP1A2, BMP2, CAMK2A, CDK5R2, CPLX1, FGF11, GALR1, GAS1, GLI3, KCNB2, KCNK3, LIME1, LRRTM1, MAPK8IP2, MLC1, PCSK9, SH3GL3, SLC26A6, SPAG5, TERT"
CESC,"ACTL7A, C9orf24, CETN2, FNDC3A, HOATZ, IGF2R, M1AP, NR0B1, PAEP, PLA2G3, PSMA8, SEPTIN7, SERPINA5, SHCBP1L, SLC2A8, SPAG8, TPGS1"
SARC,"ADIPOR1, ATAD3A, BARHL2, BDKRB1, CPNE6, CTDP1, CXCL12, CYFIP2, DRAXIN, DVL1, ENO1, FBLN5, IGFBP7, INO80, INS, KDM2B, MEAF6, NKX6_1, PLXNA1, PPARD, PSMD10, PTK2B, RGS4, SLC25A33, TMEM97, TNR, TSPYL2"
READ,"ABCC1, ADCY3, AKAP6, AQP5, AQP7, ATAT1, ATP4B, DCC, DMD, FCHO1, FERMT1, FNBP1, GP2, KCNA3, KCNB1, LRRTM3, SLC12A2, SLC22A5, SYP"
TARGET-WT,"CBLB, CDC20, DAB2, ELFN2, ERRFI1, FOXA2, HEY1, LYN, PKIA, PLN, PPP4R4, RNF34, SERPINH1, SHB, SORL1, SPRY1, TTC36"
CHOL,"JUN, PAX9, POMT2, RFPL4A, SOX15"
UCS,"CLDN8, H1_2, H1_3, H2AC6, H3C1, H3C7, KRT15, KRT83, MAP1B, MAPK8IP2, MRPL41, MRPL43, MUC17, NUP214, PKP1, TGM3"
THCA,"ALX1, BARX2, BHLHE22, E2F1, EOMES, FOXI2, GLI1, HEY2, HOXC10, HOXC13, IRX5, LMX1B, NKX2_5, NOTO, SALL3, SALL4, SNAI1, SNCA, SOX15, VAX2, ZBTB21, ZBTB8B, ZNF560, ZNF730, ZSCAN4"
DLBC,"FAM169A, FAM209A, GPER1, ITPRIP, NDC1, NOS1AP, NUCB2, PLCB1, RETSAT, RRP12, SEH1L, TMEM170A"
PAAD,"ADD3, APC2, AVIL, CAPZA1, CAPZA2, CAPZA3, CARMIL1, CARMIL2, CIB1, CLASP2, EPS8, F2RL1, GSPT1, HDAC6, KATNB1, LIMA1, MAP1B, MAP1S, MTPN, OGFOD1, PDXP, PIK3CA, SMCR8, SPTBN2, SPTBN4, TECPR1, TMEM39A, TMOD1, TMOD2, TMOD3, TPX2, TRIM54, VILL, WASHC2C, WDR1"
TARGET-AML,"ACTN1, CLASP1, CLDN19, CX3CL1, DGKE, DGKH, DST, EPHB2, F3, FOXA2, GGCX, GNAS, GPX1, ITGB3, KANK1, LNPK, MSX2, NOTCH4, PABPC4, PLCG2, PROCR, RAP2B, RHOA, SRSF6, ST3GAL4, TMEFF2, TYRO3, VTN, WFDC1, WNT7A, YAP1"
LAML,"ABCC8, ATP13A2, CALR, CARF, CPNE6, DAXX, DNMT3A, EIF2AK3, G6PD, HSD17B1, ITPKA, KCNMB1, MCOLN1, MNAT1, NEUROD2, OTC, PARP1, PEF1, PPIF, RYR1, SEC31A, SHH, SLC11A2, SLC30A10, SYT17, TCIRG1, THBS1"
KIRP,"AK3, AS3MT, CRABP2, ELOVL2, FOXA2, GADL1, GATM, GSS, INPP5J, MYH6, NPC1L1, NR5A2, PDE8B, PLAAT3, PRKN, PSAT1, PTGIS, PYCR1, RDH10, RRM2, SERPINA12, SLC52A3, SLC7A11, SULT1E1, TRIB3, TTC36, UAP1"
KIRC,"AQP1, ARHGEF1, BMP1, CHAC1, CLEC2D, EPHB2, FLRT3, FRK, GIP, GNRH1, GRB10, HHLA2, INHBE, ISG15, ITGA2, KL, MDK, NPNT, OPRD1, PDGFD, PIDD1, PLCL1, PLG, PTGER1, SEMA3E, SEMA3G, SHC1, TFF1, TLN2, TNFSF14, UCN, WNT10B, WNT9B"
TARGET-NBL,"AFDN, ATF2, CALCOCO1, CALR, CLEC5A, CPEB3, EIF4G1, EXOSC3, FOXD1, HMCES, HNRNPD, KLRC4_KLRK1, MDK, NAT8B, NSUN5, NTSR1, PIK3R1, SLC38A2, STOML2, TFRC, UCN, ZFPM1"
KICH,"SLC16A1, SLC19A2, SLC27A1, SLC38A7, SLC6A11, SLC7A10, SLC7A11, SLC7A3, SLCO3A1, SLCO6A1"
TARGET-OS,"ARNTL, CENPT, CHD1L, CREBZF, DDX21, DLX1, DLX2, EHMT2, ERCC4, FANCC, GSX2, HOXA6, HSFY1, LBX2, MAFF, MEF2A, MEF2C, MUC1, MXI1, MYC, NHEJ1, PPARG, RBM34, RFX3, UHRF2"
TARGET-ALL,"ABCB1, ACVR1B, ADGRA2, AOC3, APOA4, CALR, CD1C, CD48, CD53, CLSTN3, CYP2W1, DCBLD2, DSG2, EPHA4, ERP29, FGFR2, H1_1, HLA_DMA, IL12RB2, IL13RA1, LMAN2, MICB, MPL, NTRK1, PAM, SLC39A6, SLC9A3, TNFRSF18, TSPAN33, VAMP5"
THYM,"CRY1, CTCF, DLX2, EOMES, ESRRB, HESX1, NFIC, PBX2, PBX3, SMAD2, SMAD4, SOX5, TCF3, ZNF24, ZNF266, ZNF439, ZNF445, ZNF583, ZNF704, ZNF71, ZNF813, ZSCAN5A"
MESO,"CADPS2, CALM3, CHRNA5, CSPG5, DOC2A, DVL1, ERC2, FMR1, HTR1B, PREPL, RAP1B, SCRIB, SNAP29, STX11, STX1A, STXBP5, SYT4, SYT7, VAMP2"
LGG,"ALG6, ANG, B3GNT5, BANK1, CDC123, CHCHD1, CHST15, CPEB3, CTIF, DPY19L1, EIF3L, EIF4ENIF1, FUT11, FUT9, GALNT14, GALNT3, GALNT5, HS2ST1, HS3ST3B1, HSPB1, IGF2BP2, IGF2BP3, LARP4B, LDB1, MRPL43, MRPS16, NDST4, NOLC1, OGA, PGAP4, PIGT, PLCB1, PLOD1, PLOD3, POGLUT3, PSTK, RPL15, RPL4, RPN2, RPS24, SLC2A10, TENT5B, UGCG, ZDHHC22"
PCPG,"SLC10A4, SLC16A13, SLC17A4, SLC22A1, SLC24A4, SLC25A18, SLC2A10, SLC4A11, SLC4A4, SLC4A5, SLC4A9, SLC6A1, SLC6A13, SLC6A2, SLC6A3"
UVM,"BIRC3, CAST, CST7, DPEP1, FNIP1, GAS6, IFI16, IGBP1, KLF4, MAP2K5, MICAL1, MMP9, NLE1, PAK2, PARK7, PRDX3, PTGS2, RAF1, SERPINB9, SIAH2, THBS1, TNFAIP8, TNFSF14"
ACC,"ABCC5, AHR, AOX1, CES2, CRYZ, CYP1A2, CYP26B1, CYP2C18, CYP2D7, CYP2U1, CYP46A1, E2F1, EFTUD2, EPHX1, FMO5, GSTA1, GUK1, NAT1, NQO1, RORC, TLR3"
```



