library(shiny)
library(shinyWidgets)
library(readr)
library(DT)



Table_S89_Druggability_analysis_of_gastric_cancer_genes <- read_csv("data/druggability/Table S89. Druggability analysis of gastric cancer genes.csv")
Table_S90_Binding_effect_analysis_on_wild_and_mutant_structure_of_protein <- read_csv("data/druggability/Table S90. Binding effect analysis on wild and mutant structure of protein.csv")
sequence_func <- read_tsv("data/Functional_Study/Text S1. Protein Sequence of all cancer genes.txt")
Table_S87_List_of_domains_of_cancer_Gene <- read_csv("data/Functional_Study/domain/Table S87. List of domains of cancer Gene.csv")
Free_energy_table_from_yasara <- read_csv("data/Functional_Study/Free energy table from yasara.csv")
sequence_explore <- read_tsv("data/Exploration/Text S2. Gene sequence including selected nsSNP.txt")

Table_S2_Genes_interaction_among_the_Diseases <- read_csv("data/responsible/network/Table S2.Genes interaction among the Diseases.csv")




    

ui <- navbarPage("GasCan-2.0",
                 theme= shinythemes::shinytheme("cyborg"),
                 
                 
                 
                 
                 
                 
                 tabPanel("Home" , icon = icon("home"),
                          div(img(style="width:1200px;margin-right:15px;float: center;display: block;
  margin-left: auto;
  margin-right:auto;border: 5px solid teal;
                 
                 position: relative;
  max-width: 1310px;
  height: 180px;
  
                 
                 
  
  
  color: #f1f1f1;
  
  ",
                                  
                                  
                                  src="Cover.png")),
                          br(),
                          br(),
                          br(),
                          
                          
                            
                          div(img(style="width:580px;height:600px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right:15px;border: 5px solid #555;",
                                  
                                  src="home/fig 1.png", height = 350, width = 700)),
                          
                          
                          
                          
                          
                          includeHTML("dialog_read.html"),
                          
                          
                          
                          
                          
                          
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          
                          
                         includeHTML("home.html")
                                    
                                    
                                    
                                    
                                    
                          
                          
                          
                          
                          
                          
                          
                          
    
),
tabPanel("Responsible Genes", icon = icon("dna"),
         div(img(style="width:680px;height:520px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right:15px;border: 5px solid teal;",
                 
                 src="responsible/Cover.PNG", height = 450, width = 1310)),
         
         h4(align="center", strong("Responsible Genes")),
         br(),
         p(align="left", style = "color: #ffffff;", "To extract the most common genes involved in gastric cancer, at first, we did a comprehensive exploration of genes associated with gastric carcinogenesis. Initially, we selected a total of 40 genes implicated in gastric cancer on the basis of highest number of literature resources in the ", a(href="http://www.cancer-genetics.org/","Cancer Genetics Web database." ), "We have additionally retrieved the genes related to other types of cancer such as prostate, colon, bladder, breast and lung cancer. The interactions between GC and other cancer related genes were visualized with Cytoscape network analysis tool. Subsequently, only the genes found to be common between prostate, colon, bladder, breast, lung cancer and GC were selected to probe their interactions with one other. Following this, we utilized the GeneMania tool to reveal the representing group and their interactions with the other genes.", br(), br(), 

"We have found a group of 8 gastric cancer genes (CASP3, CD44, VEGFA, MUC1, CDKN1B, KIT, PIK3CA, TP53) which are common with bladder cancer genes. These genes share a common co-expression and pathway network.  Four genes namely, TP53, STK11, CDH1, and PTEN of GC were found to be linked with breast cancer as well. These four genes interacted with the breast cancer genes in a manner of co-expression, genetic interactions and pathway network. The interaction network revealed that TP53, PIK3CA, APC, MSH2, KRAS genes are common between GC and colon cancer. These genes were observed to interact in a way of co-expression, genetic interactions and co-localization network. We have identified only 2 genes (KRAS and MET) of GC which were found to be associated with lung cancer too. These two GC genes again showed interactions with lung cancer related genes in terms of co-expression, genetic interactions and pathway network. Interestingly, we have observed that there were no associations between GC genes with prostate cancer genes. By cross-referencing with all of the interactions data, we have identified 16 GC associated genes (CASP3, CD44, VEGFA, MUC1, CDKN1B, KIT, PIK3CA, TP53, STK11, CDH1, PTEN, PIK3CA, APC, MSH2, KRAS, KRAS, MET) that interacted with the genes related to  other types of cancer.
"),
         
         br(),
         br(),
         br(),
         br(),br(),
br(),
br(),br(),
         
         
         
         
         
         
       wellPanel( 
         
    tabsetPanel(
      
      
        tabPanel("Cancer Gene Interaction",
                 br(),
                 h3(align="center","Cancer Gene Interaction Network"),
                 tags$hr(),
                 br(),
                 br(),
                 
                 
                   
                   
                   
                     awesomeRadio("cancer_interaction","Choose Cancer Type",
                              c("Gastric & Bladder Cancer", "Gastric & Breast Cancer", "Gastric & Colon Cancer", "Gastric & Lung Cancer"), inline = TRUE, 
                              checkbox = TRUE
                              ),
                   
                 
                br(),
                br(),
                
                   
                   
                sidebarLayout(    
                       
               sidebarPanel( tableOutput("table_cancer_type")
               ),
                
                mainPanel( br(),
                           br(), 
                       imageOutput("img_cancer_type")
                
                       )
                
                     
                    
                       
                )
                       
                     
                     
                     
                       
                     
                   
                   
                   
                   
                   
                 
                 
                 
                 
                 
                 
                 
            
        ),    
        tabPanel(
            "Common Genes Among Cancer Types",
            br(),
            DTOutput("common_genes")
        )
        
        
         )
    
    
)
    
         ),
tabPanel("Exploration of nsSNP", icon = icon("microscope"),
         includeHTML("explore.html"),
         
         br(),
         tags$hr(),
         br(),
         
         div(img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
                 
                 src="Cover_explore.png", height = 180, width = 180),
         
         
        
         p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", "It is crucial to uncover whether an nsSNP could affect the protein function and contribute to the disease. The present study is based on identifying nsSNPs of GC-linked genes which might have the severely deleterious effect on its gene product. Initially, we have sorted out 44,000 nsSNPs from 1,45,342 SNPs of 40 genes of GC reported in available databases. We have then employed 7 tools namely, SIFT, PolyPhen2, PMut, MutPred, SNAP2, SNP&GO and PANTHER for the selection of nsSNPs with the most deleterious/damaging impact on the respective GC gene. At the end of this stringent filtering pipeline, only one nsSNP was chosen as the most damaging nsSNP for each gene and assumed to be involved in alteration of protein function. In total, we identified 11363 missense nsSNPs located within the 40 GC genes, 474 of whom were predicted to be damaging and 40 to be the most damaging.", br(),
           
           br(),
           

"The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. But we were unable to find out the specific restriction enzyme of each of the 40 GC genes for the specific location of polymorphism. 
")),
         br(),
         br(),
         
 wellPanel(        
         
         tabsetPanel(
           
           
           
           
           
           
           tabPanel("Primer Design",
                    
                    br(),
                    
                    
                    
                    h3(align="center", strong("RFLP & Allele Specific Primer")),
                    br(),
                   
                    p(style="color: #f1f1f1;","For the selected nsSNP association with the cancer, we have designed primer for both allele
specific PCR and PCR-RFLP assay. Allele specific PCR is effective for SNP genotyping and
mutation detection while PCR-RFLP assay is convenient for detecting any single nucleotide base
change at a specific restriction site. Firstly, we extracted respective gene sequence with the
damaging nsSNP from ", a(href ="https://ncbi.nlm.nih.gov", strong("NCBI") ), "For the designing of allele specific primer we used ", a(href="http://www.bioinformatics.nl/cgi-bin/primer3plus/primer3plus.cgi", strong("primer3plus")), "
and ", a(href="http://www.premierbiosoft.com/netprimer/", strong("NetPrimer")), "web tools Primer for PCR-RFLP assay was designed after selecting restriction enzyme that is able to cut at a specific mutation site of the gene sequence. Two web-based resources ", a(
href="http://www.labtools.us/nebcutter-v2-0/", strong("NEBcutter V2.0"))," and ", a(href="http://rna.lundberg.gu.se/cutter2/", strong("Webcutter "))," were
utilized to select restriction enzyme for PCR-RFLP. Both types of primer were designed for wild
and mutant gene sequence."),
                    
                    br(),
                    br(),  
                    
                    tabsetPanel(
                      
                      tabPanel("Allele Specific Primer",
                               
                               
                               br(),
                               sidebarLayout(
                                 
                                 sidebarPanel(  
                                   
                                   selectInput("allele_gene", "Explore Genes",
                                               c("CEACAM5","APC","CD44","CDH1","ABCB1","CTNNA1","CTNNB1","DCC","EPCAM","KIT","KITLG","CDKN1B", "BAX","BMPR1A","CASP3")
                                   )
                                 ),
                                 
                                 mainPanel(
                                   
                                   uiOutput("allele_html")
                                   
                                   
                                 )
                                 
                                 
                                 
                                 
                                 
                               )
                               
                               
                               
                               ),
                      tabPanel("PCR-RFLP Primer",
                               br(),
                               sidebarLayout(
                                 
                                 sidebarPanel(  
                                   
                                   selectInput("rflp_gene", "Explore Genes",
                                               c("ABCB1","APC","CDKN1B","CTNNB1")
                                   )
                                 ),
                                 
                                 mainPanel(
                                   
                                   uiOutput("rflp_html")
                                   
                                   
                                 )
                                 
                                 
                                 
                                 
                                 
                               )
                               
                      )
                      
                    )
             
             
             
           ),



tabPanel("In-vitro Association",
         br(),
         br(),
         br(),
         
         div(img(style="width:800px;height:450px;margin-right:15px;float: center;display: block;
  margin-left: auto;
  margin-right:auto;border: 5px solid #555;
                 box-sizing: border-box;
                 position: relative;
  max-width: 1310px;
  margin: 0 auto;
                 vertical-align: middle;
                 
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 1310;
  padding: 20px;",
                 
                 
                 src="invitro.png", height = 450, width = 1310),
         br(),br(),br(),
        div(style="background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;",
            br(),
            includeHTML("data/Exploration/Method.html")    
         
         
         
        )
         
         
         )
         
         
         
         
         
         
         
          
         
         
),




    tabPanel("Damage Prediction",
         
         br(),
         br(),
         br(),
         
         
         sidebarLayout(
           
           
           sidebarPanel(
             selectInput("html_gene", "Explore Genes",
                         c("APC","ABCB1","AURKA","BAX", "BMPR1A","CASP3","CD44","CDH1","CDKN1B","CEACAM5","CTNNA1","CTNNB1","DCC","EPCAM","FOS","KIT","KITLG","KRAS","KRT20","MALT1","MET","MGMT","MMP2","MSH2","MTHFR","MUC1","MYC","PCNA","PIK3CA","PTEN","PTGS2","SDHA","SDHB","SDHD","SMAD4","STK11","TNF","TP53","VEGFA")
             )
             
             
           ),
           mainPanel(
             
             wellPanel(style="text-align: center; 
  ",
                       
                       
                       uiOutput("html_damage")
                       
                       
             )
             
             
           )
         ),
         
         
         br(),
         br(),
         br()
         
         
         
),
           
           
           
           
           
           
           tabPanel(
             
             
             
             "Sequences Download", icon = icon("download"),
             
             br(),
             br(),
             br(),
             downloadButton(outputId = "download_nsSNP", label = "Sequences Including nsSNPs"),br(),br(),br(),br(),br(),br(),br(),br(),br()
             
             
           ),
           
           
           
           
           
           
           
           
           
         
         
         
         
         br(),
         br(),
         br(),
         br()
             
         )
)
         ),



tabPanel("Functional Study", icon = icon("atom"),
         div(img(style="width:1100px;height:450px;margin-right:15px;float: center;display: block;
  margin-left: auto;
  margin-right:auto;border: 5px solid #555;
                 box-sizing: border-box;
                 position: relative;
  max-width: 1310px;
  margin: 0 auto;
                 vertical-align: middle;
                 
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 1310;
  padding: 20px;",
           
           
           src="functional.png", height = 450, width = 1310)),
         br(),br(),br(),
        h3(style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", align="center",  strong("Functional Study")),
         br(),
         
         
                
             
               
               wellPanel(
               
               
               tabsetPanel(
                 tabPanel("3D Modelling and Free Energy",
                          br(),
                          br(),
                          br(),
                          h4(align="center", strong("3D Modelling and Free Energy Deviation Calculation")),
                          br(),
                          p(align="left", style = "color: #ffffff;", "We performed BLAST against the Protein Database (PDB) to find out the structure of the closest related proteins. We selected the closest template for each protein sequence of corresponding 40 genes. Then we built the 3D model and checked the quality of each model. Further, we replaced the amino acid (nsSNP) from wild type protein sequence and assessed the quality of the model. YASARA view mutation tool carried out the mutations (G59S, P62S, L184S, L224P, A276V, L361P, R592H, T595M, and I673T) separately of each final product of corresponding genes and it showed the decreased free energy for all the mutant models than the wild type models. These results point toward noteworthy change in the structure of each protein that can demolish its natural function."),
                          
                          br(),
                          br(),
                          tabsetPanel(
                            tabPanel("YASARA View", icon = icon("atom"),
                                     br(),
                                     br(),
                                     
                                     
                                     
                                     
                                     
                                     pickerInput("radio_func_yasara", "Select Gene From Drop Down Menu",
                                                 c("ABCB1","APC","AURKA","BAX", "BMPR1A","CASP3","CD44","CDH1","CDKN1B","CEACAM5","CTNNA1","CTNNB1","DCC","EPCAM","FOS","KIT","KITLG","KRAS","KRT20","MALT1","MET","MGMT","MMP2","MSH2","MTHFR","MUC1","MYC","PCNA","PIK3CA","PTEN","PTGS2","RUNX3","SDHA","SDHB","SDHD","SMAD4","STK11","TNF","TP53","VEGFA"),
                                                 options = list(
                                                   `live-search` = TRUE),
                                                 width = "50%",
                                                 inline = F
                                                 
                                                 
                                                 
                                     ),br(),
                                     br(),
                                     
                                     
                                     
                                     
                                     
                                     uiOutput("image_func_yasara")  , br(),br(),br(),br(),
                                     br(),
                                     br(),
                                     column(6, align="center", br(),
                                            uiOutput("image_func_yasara_download")
                                            
                                     ), 
                                     
                                     
                                     column(6, align="center", br(),
                                            uiOutput("image_func_yasara_download2")
                                            
                                     ),
                                     
                                     br(),br(),br(),br(),br(),br(),br(),br(),br(),  
                                     
                                     
                                     br(),
                                     br(),br(),
                                     br(),br(),
                                     br(),br(),
                                     br(),br(),
                                     br(),br(),
                                     br()
                                     
                                     
                                     
                                     
                                     
                                     
                                     
                            ),
                            
                            
                            tabPanel("3D Structure", icon = icon("disease"),
                                     
                                     br(),br(),br(),
                                     
                                     sidebarLayout(
                                       sidebarPanel(
                                         selectInput("radio_func_3d", "Choose Gene",
                                                      c("ABCB1","APC","AURKA","BAX", "BMPR1A","CASP3","CD44","CDH1","CDKN1B","CEACAM5","CTNNA1","CTNNB1","DCC","EPCAM","FOS","KIT","KITLG","KRAS","KRT20","MALT1","MET","MGMT","MMP2","MSH2","MTHFR","MUC1","MYC","PCNA","PIK3CA","PTEN","PTGS2(1)","PTGS2(2)","SDHA","SDHB","SDHD","SMAD4","STK11","TNF","TP53","VEGFA")
                                                      )
                                       ),
                                       
                                       mainPanel(
                                                           
                                                           imageOutput("image_func_3d")  , br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br()  
                                                           
                                       
                                       )
                                       
                                     )
                                     
                                     
                                     
                                     
                                     )
                            
                            
                            
                            
                          )
                          
                          
                          ),
                 
                 
                 tabPanel("Domain Association",
                          
                          setBackgroundImage(
                            src = "web-bg/bg.jpg"
                          ),
                          
                          
                          br(),
                          br(),
                          br(),
                          h4(align="center", strong("Domain Mapping of the nsSNPs")),
                          br(),
                          p(align="left", style = "color: #ffffff;", "To assess the impact on final protein product of 40 GC genes, we retrieved the amino acid sequences of final product (protein) of the corresponding genes. The SNPs in domain regions have been thought to be strong candidates that alter protein functions. Thus, we have tried to explore whether our identified nsSNPs fall within the domain regions. We have utilized ScanProsite, Pfam and InterPro tools to find out all possible domains of the 40 GC gene products. Our findings suggest that most of the selected nsSNP were located within some domains or motif regions. We have addressed the individual domain IDs, name of the domains, their functions as well as their specific positions within the protein."),
                          
                          br(),
                          br(),
                          tabsetPanel(
                            
                            
                            tabPanel("SNP in Domain Region of Cancer Gene", icon = icon("images"),
                                     
                                     
                                     wellPanel(style = "background-color: #ffffff;",
                                               
                                               br(),br(),
                                               img(align="center", src = "Functional/Domain/Slide1.PNG"),
                                               img(align="center", src = "Functional/Domain/Slide2.PNG"),
                                               img(align="center", src = "Functional/Domain/Slide3.PNG"),
                                               img(align="center", src = "Functional/Domain/Slide4.PNG")
                                               
                                     )
                                     
                                     
                            ),
                          
                          
                            
                            tabPanel(
                              "List of Domains of Cancer Gene",
                              
                              tags$style(HTML("
                    .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing, .dataTables_wrapper .dataTables_paginate {
                    color: #ffffff;
                    }

                    thead {
                    color: #ffffff;
                    }

                     tbody {
                    color: #000000;
                    }

                   "
                                              
                                              
                              )),
                                     br(),br(),
                              
                              DT::DTOutput("table_domain"),
                              br()
                             
                            )
                            
                            
                            
                            
                            
                            
                          )
                          ),
                 tabPanel("Free Energy Table From YASARA",
                          br(),
                          DTOutput("func_all_table"),
                          br()
                          
                          ),
                 tabPanel("Download Sequences", icon = icon("download"),
                          
                         
                          br(),
                          downloadButton(outputId = "download_func_seq", label = "Protein Sequence of All Cancer Genes"),br(),br(),br(),br(),br(),br(),br(),br(),br()
                         
                          
                          
                          )
               
               ) 
               
               
                 
             )
                 
                 
             
             
         
         ),



tabPanel("Druggability",icon = icon("pills"),
         div(img(style="width:680px;height:400px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right:15px;border: 5px solid #555;",
           
           src="druggability_illustration.jpg", height = 450, width = 1310)),
         br(),
         br(),
         h4(align="center", strong("Putative Effect of nsSNPs on Drug binding")),
         br(),
         p(align="left", style = "color: #ffffff;", "Drug binding analysis was also carried out to confirm the structural variation and possible dysfunction of the final product of each 40 GC genes between wild type and mutant model. We have utilized DrugBank to select the drugs against the protein receptor of 40 GC genes. DrugBank suggested that some drugs were available against the corresponding protein of 10 GC genes. Thereafter, we have performed molecular docking analysis between the suggested drugs and protein receptor of 10 GC genes. We have found the result of 8 GC genes which bound the different interacting residues with different binding affinity where the same docking area was used for docking runs. These results confirmed the structural variation and drug could not be effective against mutant model if the individuals with this polymorphism."),
         
         br(),
         br(),
         br(),
         br(),br(),
         br(),br(),
         br(),br(),
         br(),
         
         
             
                wellPanel(
             
                tabsetPanel(
                  
                  
                  tabPanel("Figures", icon = icon("images"),
                           
                           br(),br(),br(),
                           
                           sidebarLayout(
                             sidebarPanel(
                               awesomeRadio("radio_drug", "Choose Drug",
                                            c("Amrinone(ABCB1)","Dinoprostone","Dipyridamole", "Nelfinavir","Rifampicin","Sulfinpyrazone","Verapamil","Ethanol","Imatinib","Captopril","Fluorouracil","Methotrexate","Lenalidomide","Mafenamic-Acid","Amrinone(TNF)","Gliclazide"))
                             ),
                             
                             mainPanel(
                                                 
                                                 imageOutput("image_drug")    
                                                 
                             
                             )
                             
                           )
                  ),
                  
                  
                  
                  
                  
                    tabPanel("Druggability Analysis of Gastric Cancer Genes", 
                             
                             br(),
                             br(),
                             br(),
                             sidebarLayout(
                               
                               sidebarPanel(pickerInput("drugs", "Drugs (Multiples can be selected)",
                                                        choices = unique(Table_S89_Druggability_analysis_of_gastric_cancer_genes$`Drug Name`),
                                                        selected = unique(Table_S89_Druggability_analysis_of_gastric_cancer_genes$`Drug Name`),
                                                        options = list(`actions-box` = TRUE),
                                                        multiple = TRUE),
                                            downloadButton("download_drugdt1", "Download as You Filtered")
                              
                                            ),
                             mainPanel(DTOutput("table_drug1"),
                                       br()
                                       )
                             
                             
                             )
                             
                             ),
                    
                    
                    
                    tabPanel("Binding Effect Analysis on Wild and Mutant Structure of Protein" , br(),
                             br(),
                             br(),
                             sidebarLayout(
                               
                               sidebarPanel(pickerInput("drugs2", "Drugs (Multiples can be selected)",
                                                        choices = unique(Table_S90_Binding_effect_analysis_on_wild_and_mutant_structure_of_protein$Drug),
                                                        selected = unique(Table_S90_Binding_effect_analysis_on_wild_and_mutant_structure_of_protein$Drug),
                                                        options = list(`actions-box` = TRUE),
                                                        multiple = TRUE),
                                            downloadButton("download_drugdt2", "Download Data")
                                
                                            ),
                               mainPanel(
                                 
                                 DTOutput("table_drug2"),
                                 
                                 br()
                                    
                                         )
                               
                             )
                             )
                    
                
                )   
             
)
         ),
tabPanel("Utilities",icon = icon("external-link-alt"),
         div(
         
           
           
           
           
           
           
           
           
         img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
                 
                 src="utilities/cancer.PNG", height = 180, width = 180),
             
             
             
             p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", br(), a(href="http://www.cancerindex.org/", target = "_blank","Cancer Genetics Web", icon("external-link-alt")), br(),
               
               br(),
               
               
               "The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. For this reason, we also designed the allele specific primer by which PCR could confirm the polymorphisms of genes.
", br())),
         br(),
         
         br(),
         br(),
         
         div(
           
           
           
           
           
           
           
           
           
           img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
               
               src="utilities/ClinVar.png", height = 180, width = 180),
           
           
           
           p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", a(href="https://www.ncbi.nlm.nih.gov/clinvar/", target="_blank","ClinVar, NCBI", icon("external-link-alt")), br(),
             
             br(),
             
             
             "The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. But we were unable to find out the specific restriction enzyme of each of the 40 GC genes for the specific location of polymorphism.
")),
        br(), 
        br(),
        br(),
         div(
           
           
           
           
           
           
           
           
           
           img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
               
               src="utilities/dbsnp.png", height = 180, width = 180),
           
           
           
           p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", a(href="https://www.ncbi.nlm.nih.gov/snp/", target = "_blank", "dbSNP", icon("external-link-alt")), br(),
             
             br(),
             
             
             "The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. But we were unable to find out the specific restriction enzyme of each of the 40 GC genes for the specific location of polymorphism.
")),
        
        
        br(),
        br(),
        br(),
        
        div(
          
          
          
          
          
          
          
          
          
          img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
              
              src="utilities/embl.jpg", height = 180, width = 180),
          
          
          
          p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", a(href="https://www.ebi.ac.uk/", target = "_blank","European Bioinformatics Institute", icon("external-link-alt")), br(),
            
            br(),
            
            
            "The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. But we were unable to find out the specific restriction enzyme of each of the 40 GC genes for the specific location of polymorphism.
")),
        
        br(),
        
        br(),
        br(),
        div(
          
          
          
          
          
          
          
          
          
          img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
              
              src="utilities/hgmd.jpg", height = 180, width = 180),
          
          
          
          p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", a(href="http://www.hgmd.cf.ac.uk/ac/index.php", target = "_blank","Human Gene Mutation Database (HGMD)", icon("external-link-alt")), br(),
            
            br(),
            
            
            "The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. But we were unable to find out the specific restriction enzyme of each of the 40 GC genes for the specific location of polymorphism.
")),
        br(),
        br(),
        br(),
        div(
          
          
          
          
          
          
          
          
          
          img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
              
              src="utilities/ncbi.png", height = 180, width = 180),
          
          
          
          p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", a(href="https://www.ncbi.nlm.nih.gov/", target = "_blank","National Center for Biotechnology Information", icon("external-link-alt")), br(),
            
            br(),
            
            
            "The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. But we were unable to find out the specific restriction enzyme of each of the 40 GC genes for the specific location of polymorphism.
")),
        
        br(),
        br(),
        br(),
        div(
          
          
          
          
          
          
          
          
          img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
              
              src="utilities/omim.png", height = 180, width = 180),
          
          
          
          p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", a(href="https://www.ncbi.nlm.nih.gov/omim", target = "_blank","OMIM", icon("external-link-alt")), br(),
            
            br(),
            
            
            "The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. But we were unable to find out the specific restriction enzyme of each of the 40 GC genes for the specific location of polymorphism.
")),
        br(),
        br(),
        br(),
        
        div(
          
          
          
          
          
          
          
          
          img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
              
              src="utilities/pubmed.png", height = 180, width = 180),
          
          
          
          p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", a(href="https://pubmed.ncbi.nlm.nih.gov/", target = "_blank","PubMed", icon("external-link-alt")), 
            br(),
            
            br(),
            
            
            "The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. But we were unable to find out the specific restriction enzyme of each of the 40 GC genes for the specific location of polymorphism. 
")),
        
        br(),
        br(),
        br(),
        
        
        div(
          
          
          
          
          
          
          
          
          
          img(style="width:280px;height:180px;margin-right:15px;float: left;display: block;
  margin-left: 15px;
  margin-right:15px;
  
  margin-top: 15px;
  margin-bottom: 5px;
                 ",
              
              src="utilities/sib.jpg", height = 180, width = 180),
          
          
          
          p(align="left", style="
  bottom: 0;
  background: rgb(0, 0, 0); /* Fallback color */
  background: rgba(0, 0, 0, 0.5); /* Black background with 0.5 opacity */
  color: #f1f1f1;
  width: 100%;
  padding: 20px;", a(href="https://www.sib.swiss/",target="_blank","Swiss Institute of Bioinformatics", icon("external-link-alt")), br(),
            
            br(),
            
            
            "The current study also set out to further a comprehensive framework for robust in vitro detection of association of our target nsSNP with 40 GC genes. Therefore, we have predicted a primer set (forward primer, reverse primer), allele specific primer and restriction enzyme to accelerate the association study where the human blood sample have to be utilized for DNA extraction. We have retrieved the gene sequence (length 500-100bp) of each gene where the nsSNP was included within this sequence. We also checked the similarity of the designed primers against Human Genome using the BLAST tool. We have scrutinized the primer to make sure that there would be no non-specific amplification. To identify the specific polymorphism (nsSNP) within the respective gene sequence from PCR product we have also selected the most suitable restriction enzyme. But we were unable to find out the specific restriction enzyme of each of the 40 GC genes for the specific location of polymorphism.
")),
         
         
         
         
         
         
         br()
         
         
         
         ),








tabPanel("Contact", icon = icon("id-card"),
         HTML("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"),
         
         br(),
         includeHTML("test.html"),
         
         br(),br(),br(),br(),
         h4(align="center",strong("Give Us Your Feedback")),
         p(style = "color: #ffffff;","The success of GasCanBase 2.0 depends on the satisfaction of our visitors. We welcome your questions and suggestions as to how we can improve our application.

GasCanBase 2.0 is operated by the Bioinfomatics Division, under The National Institute of Biotechnology, Ministry of Science & Technology")
         
         
         )
)





server <- function(input, output) {
    
    observeEvent(input$about, {
       showModal(modalDialog(title= "About GasCanBase-2.0", dialog))})
       
    observeEvent(input$about_res, {
        showModal(modalDialog(title= "Putative Effect of nsSNPs on Drug binding", dialog2))})
    
    
    
    filtered_drug <- reactive({
       
        data_drug <- Table_S89_Druggability_analysis_of_gastric_cancer_genes
        data_drug <- subset(
            data_drug,
            `Drug Name` %in% input$drugs)
        
          
    })
    
    
    filtered_data <- reactive({
        
        data <- gapminder
        data <- subset(
            data,
            continent %in% input$continents &
                year >= input$years[1] & year <= input$years[2] & lifeExp >= input$life[1] & lifeExp <= input$life[2])
        
    })
    
    output$download_data <- downloadHandler(
        filename = "Emon's Responsive Compliment.csv",
        content = function(file) {
           
            data <- filtered_data()
            write.csv(data, file, row.names = FALSE)
        }
    )
    
    
   output$common_genes <- renderDT({
     data <- Table_S2_Genes_interaction_among_the_Diseases
     data
     
   }) 
    
    
    
    
    
    
    output$download_func_seq <- downloadHandler(
      filename = "Protein Sequence of All Cancer Genes.txt",
      content = function(file) {
        
        data <- sequence_func
        write_tsv(data, file)
      }
    )
    
    
    output$download_nsSNP <- downloadHandler(
      filename = "Sequences Including nsSNPs.txt",
      content = function(file) {
        
        data <- sequence_explore
        write_tsv(data, file)
      }
    )
    
    
    
    
    output$image_drug <- renderImage({
      
       
      
      
      if(input$radio_drug == "Amrinone(ABCB1)"){            
        list(src = "www/figures_drugs/1.PNG",height = 400, width = 800)
      }                                        
      else if(input$radio_drug == "Dinoprostone"){
        list(src = "www/figures_drugs/2.PNG", height = 400, width = 800)
      }
      else if(input$radio_drug == "Dipyridamole"){
        list(src = "www/figures_drugs/3.PNG", height = 400, width = 800)
      }
      else if(input$radio_drug == "Nelfinavir"){
        list(src = "www/figures_drugs/4.PNG", height = 400, width = 800)
      }
      else if(input$radio_drug == "Rifampicin"){
        list(src = "www/figures_drugs/5.PNG", height = 400, width = 800)
      }
      else if(input$radio_drug == "Sulfinpyrazone"){
        list(src = "www/figures_drugs/6.PNG", height = 419, width = 800)
      }
      else if(input$radio_drug == "Verapamil"){
        list(src = "www/figures_drugs/7.PNG", height = 400, width = 800)
      }
      else if(input$radio_drug == "Ethanol"){
        list(src = "www/figures_drugs/8.PNG", height = 415, width = 800)
      }
      else if(input$radio_drug == "Imatinib"){
        list(src = "www/figures_drugs/9.PNG", height = 418, width = 800)
      }
      else if(input$radio_drug == "Captopril"){
        list(src = "www/figures_drugs/10.PNG", height = 400, width = 800)
      }
      else if(input$radio_drug == "Fluorouracil"){
        list(src = "www/figures_drugs/11.PNG", height = 400, width = 800)
      }
      else if(input$radio_drug == "Methotrexate"){
        list(src = "www/figures_drugs/12.PNG", height = 418, width = 800)
      }
      else if(input$radio_drug == "Lenalidomide"){
        list(src = "www/figures_drugs/13.PNG", height = 400, width = 800)
      }
      else if(input$radio_drug == "Mafenamic-Acid"){
        list(src = "www/figures_drugs/14.PNG", height = 416, width = 800)
      }
      else if(input$radio_drug == "Amrinone(TNF)"){
        list(src = "www/figures_drugs/15.PNG", height = 400, width = 800)
      }
      else if(input$radio_drug == "Gliclazide"){
        list(src = "www/figures_drugs/16.PNG", height = 400, width = 800)
      }
      
      
      
      
      
      
    })
    
    
    output$image_func_3d <- renderImage({  
      
      
      if(input$radio_func_3d == "ABCB1"){            
        list(src = "www/Functional/3d_structure/Slide1.PNG",height = 550, width = 800)
      }                                        
      else if(input$radio_func_3d == "APC"){
        list(src = "www/Functional/3d_structure/Slide2.PNG", height =550, width = 800)
      }
      else if(input$radio_func_3d == "AURKA"){
        list(src = "www/Functional/3d_structure/Slide3.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "BAX"){
        list(src = "www/Functional/3d_structure/Slide4.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "BMPR1A"){
        list(src = "www/Functional/3d_structure/Slide5.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "CASP3"){
        list(src = "www/Functional/3d_structure/Slide6.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "CD44"){
        list(src = "www/Functional/3d_structure/Slide7.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "CDH1"){
        list(src = "www/Functional/3d_structure/Slide8.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "CDKN1B"){
        list(src = "www/Functional/3d_structure/Slide9.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "CEACAM5"){
        list(src = "www/Functional/3d_structure/Slide10.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "CTNNA1"){
        list(src = "www/Functional/3d_structure/Slide11.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "CTNNB1"){
        list(src = "www/Functional/3d_structure/Slide12.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "DCC"){
        list(src = "www/Functional/3d_structure/Slide13.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "EPCAM"){
        list(src = "www/Functional/3d_structure/Slide14.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "FOS"){
        list(src = "www/Functional/3d_structure/Slide15.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "KIT"){
        list(src = "www/Functional/3d_structure/Slide16.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "KITLG"){
        list(src = "www/Functional/3d_structure/Slide17.PNG", height =550, width = 800)
      }
      else if(input$radio_func_3d == "KRAS"){
        list(src = "www/Functional/3d_structure/Slide18.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "KRT20"){
        list(src = "www/Functional/3d_structure/Slide19.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "MALT1"){
        list(src = "www/Functional/3d_structure/Slide20.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "MET"){
        list(src = "www/Functional/3d_structure/Slide21.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "MGMT"){
        list(src = "www/Functional/3d_structure/Slide22.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "MMP2"){
        list(src = "www/Functional/3d_structure/Slide23.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "MSH2"){
        list(src = "www/Functional/3d_structure/Slide24.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "MTHFR"){
        list(src = "www/Functional/3d_structure/Slide25.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "MUC1"){
        list(src = "www/Functional/3d_structure/Slide26.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "MYC"){
        list(src = "www/Functional/3d_structure/Slide27.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "PCNA"){
        list(src = "www/Functional/3d_structure/Slide28.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "PIK3CA"){
        list(src = "www/Functional/3d_structure/Slide29.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "PTEN"){
        list(src = "www/Functional/3d_structure/Slide30.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "PTGS2(1)"){
        list(src = "www/Functional/3d_structure/Slide31.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "PTGS2(2)"){
        list(src = "www/Functional/3d_structure/Slide32.PNG", height =550, width = 800)
      }
      else if(input$radio_func_3d == "SDHA"){
        list(src = "www/Functional/3d_structure/Slide33.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "SDHB"){
        list(src = "www/Functional/3d_structure/Slide34.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "SDHD"){
        list(src = "www/Functional/3d_structure/Slide35.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "SMAD4"){
        list(src = "www/Functional/3d_structure/Slide36.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "STK11"){
        list(src = "www/Functional/3d_structure/Slide37.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "TNF"){
        list(src = "www/Functional/3d_structure/Slide38.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "TP53"){
        list(src = "www/Functional/3d_structure/Slide39.PNG", height = 550, width = 800)
      }
      else if(input$radio_func_3d == "VEGFA"){
        list(src = "www/Functional/3d_structure/Slide40.PNG", height = 550, width = 800)
      }
      
      
      
    })
    
    
    
    
    
    
    output$image_func_yasara <- renderUI({  
      
      
      if(input$radio_func_yasara == "ABCB1"){            
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
          src = "Functional/Free energy calculation/ABCB1/yasara/ABCB1_W.png" 
          
          
          
          ),
                 img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
           src = "Functional/Free energy calculation/ABCB1/yasara/I736K.png"))
      }
      
      
      else if(input$radio_func_yasara == "APC"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/APC/yasara/APC-W.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/APC/yasara/D1841Y.png"))
      }
      
      
      
      else if(input$radio_func_yasara == "AURKA"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/AURKA/yasara/PIC-W.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/AURKA/yasara/PIC-G325W.png"))
      }
      else if(input$radio_func_yasara == "BAX"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/BAX/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/BAX/yasara/NEW.png"))
      }
      
      
      
      else if(input$radio_func_yasara == "BMPR1A"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/BMPR1A/yasara/PIC-W.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/BMPR1A/yasara/PIC-R443C.png"))
      }
      else if(input$radio_func_yasara == "CASP3"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/CASP3/yasara/PIC-W.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/CASP3/yasara/PIC-S249Y.png"))
      }
      else if(input$radio_func_yasara == "CD44"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/CD44/yasara/W.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/CD44/yasara/R673Q.png"))
      }
      else if(input$radio_func_yasara == "CDH1"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/CDH1/yasara/PIC-W.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/CDH1/yasara/PIC-V832M.png"))
      }
      else if(input$radio_func_yasara == "CDKN1B"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/CDKN1B/yasara/wild.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/CDKN1B/yasara/PIC-R15W.png"))
      }
      else if(input$radio_func_yasara == "CEACAM5"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/CEACAM5/yasara/wild.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/CEACAM5/yasara/Q137P.png"))
      }
      else if(input$radio_func_yasara == "CTNNA1"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/CTNNA1/yasara/wild.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/CTNNA1/yasara/R194K.png"))
      }
      else if(input$radio_func_yasara == "CTNNB1"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/CTNNB1/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/CTNNB1/yasara/G34E.png"))
      }
      else if(input$radio_func_yasara == "DCC"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/DCC/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/DCC/yasara/NEW.png"))
      }
      else if(input$radio_func_yasara == "EPCAM"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/EPCAM/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/EPCAM/yasara/I277M.png"))
      }
      else if(input$radio_func_yasara == "FOS"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/FOS/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/FOS/yasara/V77G.png"))
      }
      else if(input$radio_func_yasara == "KIT"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/KIT/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/KIT/yasara/D816V.png"))
      }
      else if(input$radio_func_yasara == "KITLG"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/KITLG/yasara/wild.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/KITLG/yasara/d149y.png"))
      }
      else if(input$radio_func_yasara == "KRAS"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/KRAS/yasara/wild.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/KRAS/yasara/G60R.png"))
      }
      else if(input$radio_func_yasara == "KRT20"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/KRT20/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/KRT20/yasara/H208Y.png"))
      }
      else if(input$radio_func_yasara == "MALT1"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/MALT1/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/MALT1/yasara/I60G.png"))
      }
      else if(input$radio_func_yasara == "MET"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/MET/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/MET/yasara/D1265Y.png"))
      }
      else if(input$radio_func_yasara == "MGMT"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/MGMT/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/MGMT/yasara/R159Q.png"))
      }
      else if(input$radio_func_yasara == "MMP2"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/MMP2/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/MMP2/yasara/C288F.png"))
      }
      else if(input$radio_func_yasara == "MSH2"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/MSH2/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/MSH2/yasara/G751R.png"))
      }
      else if(input$radio_func_yasara == "MTHFR"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/MTHFR/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/MTHFR/yasara/R377C.png"))
      }
      else if(input$radio_func_yasara == "MUC1"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/MUC1/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/MUC1/yasara/V126M.png"))
      }
      else if(input$radio_func_yasara == "MYC"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/MYC/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/MYC/yasara/P74A.png"))
      }
      else if(input$radio_func_yasara == "PCNA"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/PCNA/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/PCNA/yasara/S39R.png"))
      }
      else if(input$radio_func_yasara == "PIK3CA"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/PIK3CA/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/PIK3CA/yasara/R555K.png"))
      }
      else if(input$radio_func_yasara == "PTEN"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/PTEN/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/PTEN/yasara/D252G.png"))
      }
      else if(input$radio_func_yasara == "PTGS2"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/PTGS2/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/PTGS2/yasara/E488G.png"))
      }
      else if(input$radio_func_yasara == "RUNX3"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/RUNX3/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/RUNX3/yasara/G335R.png"))
      }
      else if(input$radio_func_yasara == "SDHA"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/SDHA/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/SDHA/yasara/R600Q.png"))
      }
      else if(input$radio_func_yasara == "SDHB"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/SDHB/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/SDHB/yasara/C101Y.png"))
      }
      else if(input$radio_func_yasara == "SDHD"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/SDHD/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/SDHD/yasara/H102L.png"))
      }
      else if(input$radio_func_yasara == "SMAD4"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/SMAD4/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/SMAD4/yasara/R361C.png"))
      }
      else if(input$radio_func_yasara == "STK11"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/STK11/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/STK11/yasara/D194N.png"))
      }
      else if(input$radio_func_yasara == "TNF"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/TNF/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/TNF/yasara/I194N.png"))
      }
      else if(input$radio_func_yasara == "TP53"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/TP53/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/TP53/yasara/R282W.png"))
      }
      else if(input$radio_func_yasara == "VEGFA"){
        tags$div(img(style="width:630px;height:460px;margin-right:15px;float: left;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
                     src = "Functional/Free energy calculation/VEGFA/yasara/WILD.png" 
                     
                     
                     
        ),
        img(style="width:630px;height:460px;margin-right:15px;float: right;display: block;
  margin-left: auto;
  margin-right: auto; background-image: linear-gradient(rgba(0, 0, 0, 0.5), rgba(0, 0, 0, 0.5)),
  background-position: center;
  background-repeat: no-repeat;
  background-size: cover;
  position: relative;border: 5px solid black;",
            src = "Functional/Free energy calculation/VEGFA/yasara/R288W.png"))
      }
      
      
      
    })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    output$download_drugdt1 <- downloadHandler(
      filename = "Druggability Analysis of Gastric Cancer Genes.csv",
      content = function(file) {
        
        data <- filtered_drug()
        write.csv(data, file, row.names = FALSE)
      }
    )
    
    
    
    
    
    
    output$img_cancer_type <- renderImage({
      
      
      
      
      if(input$cancer_interaction == "Gastric & Bladder Cancer"){            
        list(style="border: 5px solid teal;
                 box-sizing: border-box;",src = "data/responsible/Cancer Gene Interactions Network/Gene networking of Gastric and Bladder Cancer genes.png",height = 850, width = 750)
      }                                        
      else if(input$cancer_interaction == "Gastric & Breast Cancer"){
        list(style="border: 5px solid teal;
                 box-sizing: border-box;",src = "data/responsible/Cancer Gene Interactions Network/Gene networking of Gastric and Breast Cancer genes.png", height = 850, width = 750)
      }
      else if(input$cancer_interaction == "Gastric & Colon Cancer"){
        list(style="border: 5px solid teal;
                 box-sizing: border-box;",src = "data/responsible/Cancer Gene Interactions Network/Gene networking of Gastric and Colon Cancer genes.png", height = 850, width = 750)
      }
      else if(input$cancer_interaction == "Gastric & Lung Cancer"){
        list(style="border: 5px solid teal;
                 box-sizing: border-box;",src = "data/responsible/Cancer Gene Interactions Network/Gene networking of Gastric and Lung Cancer genes.png", height = 850, width = 750)
      }
      
    })
    
    
    
    
    output$table_cancer_type <- renderTable({
      
      
      
      
      if(input$cancer_interaction == "Gastric & Bladder Cancer"){            
        data <- read.csv("data/responsible/Cancer Gene Interactions Network/Gene networking of Gastric and Bladder Cancer genes.csv")
        data
      }                                        
      else if(input$cancer_interaction == "Gastric & Breast Cancer"){
        data <- read.csv("data/responsible/Cancer Gene Interactions Network/Gene networking of Gastric and Breast Cancer genes.csv")
        data
      }
      else if(input$cancer_interaction == "Gastric & Colon Cancer"){
        data <- read.csv("data/responsible/Cancer Gene Interactions Network/Gene networking of Gastric and Colon Cancer genes.csv")
        data
      }
      else if(input$cancer_interaction == "Gastric & Lung Cancer"){
        data <- read.csv("data/responsible/Cancer Gene Interactions Network/Gene networking of Gastric and Lung Cancer genes.csv")
        data
      }
      
    })
    
    
    
    
    
    output$download_drugdt2 <- downloadHandler(
      filename = "Binding Effect Analysis on Wild and Mutant Structure of Protein.csv",
      content = function(file) {
        
        data <- Table_S90_Binding_effect_analysis_on_wild_and_mutant_structure_of_protein
        write.csv(data, file, row.names = FALSE)
      }
    )
    
    
    
    
    
    
    
    
    output$rflp_html <- renderUI({  
      
      
      if(input$rflp_gene == "ABCB1"){            
        includeHTML("data/Exploration/4. Primer design/1. PCR-RFLP Primers/Table S5.1. Primer and Restriction enzyme selection of ABCB1 gene.html")
      }                                        
      else if(input$rflp_gene == "APC"){
        includeHTML("data/Exploration/4. Primer design/1. PCR-RFLP Primers/Table S5.2. Primer and Restriction enzyme selection of APC gene.html")
      }
      else if(input$rflp_gene == "CDKN1B"){
        includeHTML("data/Exploration/4. Primer design/1. PCR-RFLP Primers/Table S5.3. Primer and Restriction enzyme selection of CDKN1B gene.html")
      }
      else if(input$rflp_gene == "CTNNB1"){
        includeHTML("data/Exploration/4. Primer design/1. PCR-RFLP Primers/Table S5.4. Primer and Restriction enzyme selection of CTNNB1 gene.html")
      }
      
      
    })
    
    
    
    
    
    output$allele_html <- renderUI({  
      
      
      if(input$allele_gene == "ABCB1"){            
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.1. Allele specific primer design on selected nsSNP of ABCB1 gene.html")
      }                                        
      else if(input$allele_gene == "CEACAM5"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.9. Allele specific primer design on selected nsSNP of CEACAM5 gene.html")
      }
      else if(input$allele_gene == "CDKN1B"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.8. Allele specific primer design on selected nsSNP of CDKN1B gene.html")
      }
      else if(input$allele_gene == "KITLG"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.15. Allele specific primer design on selected nsSNP of KITLG gene.html")
      }
      else if(input$allele_gene == "KIT"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.14. Allele specific primer design on selected nsSNP of KIT gene.html")
      }
      else if(input$allele_gene == "EPCAM"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.13. Allele specific primer design on selected nsSNP of EPCAM gene.html")
      }
      else if(input$allele_gene == "DCC"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.12. Allele specific primer design on selected nsSNP of DCC gene.html")
      }
      else if(input$allele_gene == "CTNNA1"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.10. Allele specific primer design on selected nsSNP of CTNNA1 gene.html")
      }
      else if(input$allele_gene == "CDH1"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.7. Allele specific primer design on selected nsSNP of CDH1 gene.html")
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/CDH1Allele specific primer design on selected nsSNP of  gene.html")
      }
      else if(input$allele_gene == "CD44"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.6. Allele specific primer design on selected nsSNP of CD44 gene.html")
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/CD44  Allele specific primer design on selected nsSNP of gene.html")
      }
      else if(input$allele_gene == "CASP3"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.5. Allele specific primer design on selected nsSNP of CASP3 gene.html")
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/CASP3 Allele specific primer design on selected nsSNP of gene.html")
      }
      else if(input$allele_gene == "BMPR1A"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.4. Allele specific primer design on selected nsSNP of BMPR1A gene.html")
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/BMPR1A Allele specific primer design on selected nsSNP of BMPR1A gene.html")
      }
      else if(input$allele_gene == "CTNNB1"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.11. Allele specific primer design on selected nsSNP of CTNNB1 gene.html")
      }
      else if(input$allele_gene == "APC"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.2. Allele specific primer design on selected nsSNP of APC gene.html")
      }
      else if(input$allele_gene == "BAX"){
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/Table S6.3. Allele specific primer design on selected nsSNP of BAX gene.html")
        includeHTML("data/Exploration/4. Primer design/2. Allele specific Primers/BAX Allele specific primer design on selected nsSNP of BAX gene.html")
      }
       
      
      
    })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    output$html_damage <- renderUI({  
      
      
      if(input$html_gene == "ABCB1"){            
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S47. Prediction of damaging effect on ABCB1.html")
      }                                        
      else if(input$html_gene == "APC"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S48. Prediction of damaging effect on APC.html")
      }
      else if(input$html_gene == "AURKA"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S49. Prediction of damaging effect on AURKA.html")
      }
      else if(input$html_gene == "BAX"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S50. Prediction of damaging effect on BAX.html")
      }
      else if(input$html_gene == "BMPR1A"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S51. Prediction of damaging effect on BMPR1A.html")
      }
      else if(input$html_gene == "CASP3"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S52. Prediction of damaging effect on CASP3.html")
      }
      else if(input$html_gene == "CD44"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S53. Prediction of damaging effect on CD44.html")
      }
      else if(input$html_gene == "CDH1"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S54. Prediction of damaging effect on CDH1.html")
      }
      else if(input$html_gene == "CDKN1B"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S55. Prediction of damaging effect on CDK1B.html")
      }
      else if(input$html_gene == "CEACAM5"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S56. Prediction of damaging effect on CEACAM5.html")
      }
      else if(input$html_gene == "CTNNA1"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S57. Prediction of damaging effect on CTNNA1.html")
      }
      else if(input$html_gene == "CTNNB1"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S58. Prediction of damaging effect on CTNNB1.html")
      }
      else if(input$html_gene == "DCC"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S59. Prediction of damaging effect on DCC.html")
      }
      else if(input$html_gene == "EPCAM"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S60. Prediction of damaging effect on EPCAM.html")
      }
      else if(input$html_gene == "FOS"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S61. Prediction of damaging effect on FOS.html")
      }
      else if(input$html_gene == "KIT"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S62. Prediction of damaging effect on KIT.html")
      }
      else if(input$html_gene == "KITLG"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S63. Prediction of damaging effect on KITLG.html")
      }
      else if(input$html_gene == "KRAS"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S64. Prediction of damaging effect on KRAS.html")
      }
      else if(input$html_gene == "KRT20"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S65. Prediction of damaging effect on KRT20.html")
      }
      else if(input$html_gene == "MALT1"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S66.  Prediction of damaging effect on MALT1.html")
      }
      else if(input$html_gene == "MET"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S67. Prediction of damaging effect on MET.html")
      }
      else if(input$html_gene == "MGMT"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S68. Prediction of damaging effect on MGMT.html")
      }
      else if(input$html_gene== "MMP2"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S69. Prediction of damaging effect on MMP2.html")
      }
      else if(input$html_gene == "MSH2"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S70. Prediction of damaging effect on MSH2.html")
      }
      else if(input$html_gene == "MTHFR"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S71. Prediction of damaging effect on MTHFR.html")
      }
      else if(input$html_gene == "MUC1"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S72. Prediction of damaging effect on MUC1.html")
      }
      else if(input$html_gene == "MYC"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S73. Prediction of damaging effect on MYC.html")
      }
      else if(input$html_gene == "PCNA"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S74. Prediction of damaging effect on PCNA.html")
      }
      else if(input$html_gene == "PIK3CA"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S75. Prediction of damaging effect on PIK3CA.html")
      }
      else if(input$html_gene == "PTEN"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S76. Prediction of damaging effect on PTEN.html")
      }
      else if(input$html_gene == "PTGS2"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S77. Prediction of damaging effect on PTGS2.html")
      }
     
      else if(input$html_gene == "SDHA"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S79. Prediction of damaging effect on SDHA.html")
      }
      else if(input$html_gene == "SDHB"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S80. Prediction of damaging effect on SDHB.html")
      }
      else if(input$html_gene == "SDHD"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S81. Prediction of damaging effect on SDHD.html")
      }
      else if(input$html_gene == "SMAD4"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S82. Prediction of damaging effect on SMAD4.html")
      }
      else if(input$html_gene == "STK11"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S83. Prediction of damaging effect on STK11.html")
      }
      else if(input$html_gene== "TNF"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S84. Prediction of damaging effect on TNF.html")
      }
      else if(input$html_gene == "TP53"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S85. Prediction of damaging effect on TP53.html")
      }
      else if(input$html_gene == "VEGFA"){
        includeHTML("data/Exploration/3. nsSNP damage prediction/Table S86. Prediction of damaging effect on VEGFA.html")
      }
      
      
      
    })
    
    
    
    
    
    
    
    
    
    
    output$table <- DT::renderDT({
       
        data <- filtered_data()
        data
        
    })
    
    
    output$table_domain <- DT::renderDT({
      
      data <- Table_S87_List_of_domains_of_cancer_Gene
      data
      
    })
    
    
    
    
    
    
    
    
    
    
    
    output$image_func_yasara_download <- renderUI({  
      
      
      if(input$radio_func_yasara == "ABCB1"){            
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/ABCB1/yasara/00017497-W.sce", download = "00017497-W.sce", icon("download"), " Wild-Type ABCB1 .sce File")
                    
        )
            
      }
      
      
      else if(input$radio_func_yasara == "APC"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/APC/yasara/00017630-W_303.sce", download = "00017630-W_303.sce", icon("download"), " Wild-Type APC .sce File")
                    
        )
            
      }
      
      
      
      else if(input$radio_func_yasara == "AURKA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/AURKA/yasara/00016606-aurka_wild.sce", download = "00016606-aurka_wild.sce", icon("download"), " Wild-Type AURKA .sce File")
                    
        )
            
      }
      else if(input$radio_func_yasara == "BAX"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/BAX/yasara/00017505-BAX-W.sce", download = "00017505-BAX-W.sce", icon("download"), " Wild-Type BAX .sce File")
                    
        )
            
      }
      
      
      
      else if(input$radio_func_yasara == "BMPR1A"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/BMPR1A/yasara/00017513-W.sce", download = "00017513-W.sce", icon("download"), " Wild-Type BMPR1A .sce File")
                    
        )
            
      }
      else if(input$radio_func_yasara == "CASP3"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CASP3/yasara/00017514-W.sce", download = "00017514-W.sce", icon("download"), " Wild-Type CASP3 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CD44"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CD44/yasara/00017557-W.sce", download = "00017557-W.sce", icon("download"), " Wild-Type CD44 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CDH1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CDH1/yasara/00017555-W.sce", download = "00017555-W.sce", icon("download"), " Wild-Type CDH1 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CDKN1B"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CDKN1B/yasara/00017553-W.sce", download = "00017553-W.sce", icon("download"), " Wild-Type CDKN1B .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CEACAM5"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CEACAM5/yasara/00017552-W.sce", download = "00017552-W.sce", icon("download"), " Wild-Type CEACAM5 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CTNNA1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CTNNA1/yasara/00017554-W.sce", download = "00017554-W.sce", icon("download"), " Wild-Type CTNNA1 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CTNNB1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CTNNB1/yasara/00017565-W.sce", download = "00017565-W.sce", icon("download"), " Wild-Type CTNNB1 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "DCC"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/DCC/yasara/00017563-W.sce", download = "00017563-W.sce", icon("download"), " Wild-Type DCC .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "EPCAM"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/EPCAM/yasara/00017571-W.sce", download = "00017571-W.sce", icon("download"), " Wild-Type EPCAM .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "FOS"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/FOS/yasara/00016604-W.sce", download = "00016604-W.sce", icon("download"), " Wild-Type FOS .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "KIT"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/KIT/yasara/00017575-W_PROB.sce", download = "00017575-W_PROB.sce", icon("download"), " Wild-Type KIT .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "KITLG"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/KITLG/yasara/00017578-W.sce", download = "00017578-W.sce", icon("download"), " Wild-Type KITLG .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "KRAS"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/KRAS/yasara/00017580-W.sce", download = "00017580-W.sce", icon("download"), " Wild-Type KRAS .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "KRT20"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/KRT20/yasara/00016600-wild.sce", download = "00016600-wild.sce", icon("download"), " Wild-Type KRT20 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MALT1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MALT1/yasara/00017582-W.sce", download = "00017582-W.sce", icon("download"), " Wild-Type MALT1 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MET"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MET/yasara/00017594-W.sce", download = "00017594-W.sce", icon("download"), " Wild-Type MET .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MGMT"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MGMT/yasara/00017592-W.sce", download = "00017592-W.sce", icon("download"), " Wild-Type MGMT .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MMP2"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MMP2/yasara/00017590-W.sce", download = "00017590-W.sce", icon("download"), " Wild-Type MMP2 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MSH2"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MSH2/yasara/00017639-W.sce", download = "00017639-W.sce", icon("download"), " Wild-Type MSH2 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MTHFR"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MTHFR/yasara/00017597-W.sce", download = "00017597-W.sce", icon("download"), " Wild-Type MTHFR .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MUC1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MUC1/yasara/00017598-W.sce", download = "00017598-W.sce", icon("download"), " Wild-Type MUC1 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MYC"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MYC/yasara/00017605-W.sce", download = "00017605-W.sce", icon("download"), " Wild-Type MYC .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "PCNA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/PCNA/yasara/00017603-W.sce", download = "00017603-W.sce", icon("download"), " Wild-Type PCNA .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "PIK3CA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/PIK3CA/yasara/00016598-wild.sce", download = "00016598-wild.sce", icon("download"), " Wild-Type PIK3CA .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "PTEN"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/PTEN/yasara/00017637-W.sce", download = "00017637-W.sce", icon("download"), " Wild-Type PTEN .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "PTGS2"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/PTGS2/yasara/00017615-W.sce", download = "00017615-W.sce", icon("download"), " Wild-Type PTGS2 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "RUNX3"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/RUNX3/yasara/00017614-W.sce", download = "00017614-W.sce", icon("download"), " Wild-Type RUNX3 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "SDHA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/SDHA/yasara/00017612-W.sce", download = "00017612-W.sce", icon("download"), " Wild-Type SDHA .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "SDHB"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/SDHB/yasara/00017610-W.sce", download = "00017610-W.sce", icon("download"), " Wild-Type SDHB .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "SDHD"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/SDHD/yasara/00017629-w.sce", download = "00017629-w.sce", icon("download"), " Wild-Type SDHD .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "SMAD4"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/SMAD4/yasara/00017601-W.sce", download = "00017601-W.sce", icon("download"), " Wild-Type SMAD4 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "STK11"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/STK11/yasara/00017627-W.sce", download = "00017627-W.sce", icon("download"), " Wild-Type STK11 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "TNF"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/TNF/yasara/00017624-W.sce", download = "00017624-W.sce", icon("download"), " Wild-Type TNF .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "TP53"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/TP53/yasara/00017642-W.sce", download = "00017642-W.sce", icon("download"), " Wild-Type TP53 .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "VEGFA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/VEGFA/yasara/00017588-W.sce", download = "00017588-W.sce", icon("download"), " Wild-Type VEGFA .sce File")
                    
        )
      }
      
      
      
    })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
      
    
    output$table_drug1 <- renderDT({
        
        data_drug <- filtered_drug()
        data_drug
        
    }) 
    
    
    output$func_all_table <- renderDT({
      
      data_drug <- Free_energy_table_from_yasara
      data_drug
      
    }) 
    
    
    output$table_drug2 <- renderDT({
        
        data_drug2 <- Table_S90_Binding_effect_analysis_on_wild_and_mutant_structure_of_protein
        data_drug2 <- subset(
            data_drug2,
            Drug %in% input$drugs2)
        data_drug2
        
    }) 
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    output$image_func_yasara_download2 <- renderUI({  
      
      
      if(input$radio_func_yasara == "ABCB1"){            
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/ABCB1/yasara/00017498-I736K.sce", download = "00017498-I736K.sce", icon("download"), "Mutant ABCB1 (I736K) .sce File")
                    
        )
        
      }
      
      
      else if(input$radio_func_yasara == "APC"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/APC/yasara/00017640-D1841Y.sce", download = "00017640-D1841Y.sce", icon("download"), "Mutant APC (D1841Y) .sce File")
                    
        )
        
      }
      
      
      
      else if(input$radio_func_yasara == "AURKA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/AURKA/yasara/00016605_G325W.sce", download = "00016605_G325W.sce", icon("download"), "Mutant AURKA (G325W) .sce File")
                    
        )
        
      }
      else if(input$radio_func_yasara == "BAX"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/BAX/yasara/00017504-BAX-N.sce", download = "00017504-BAX-N.sce", icon("download"), "Mutant BAX .sce File")
                    
        )
        
      }
      
      
      
      else if(input$radio_func_yasara == "BMPR1A"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/BMPR1A/yasara/00017508-R443C.sce", download = "00017508-R443C.sce", icon("download"), "Mutant BMPR1A (R443C) .sce File")
                    
        )
        
      }
      else if(input$radio_func_yasara == "CASP3"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CASP3/yasara/00017515-S249Y.sce", download = "00017515-S249Y.sce", icon("download"), "Mutant CASP3 (S249Y) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CD44"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CD44/yasara/00017556-R673Q.sce", download = "00017556-R673Q.sce", icon("download"), "Mutant CD44 (R673Q) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CDH1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CDH1/yasara/00017569-V832M.sce", download = "00017569-V832M.sce", icon("download"), "Mutant CDH1 (V832M) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CDKN1B"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CDKN1B/yasara/00017568-R15W.sce", download = "00017568-R15W.sce", icon("download"), "Mutant CDKN1B (R15W) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CEACAM5"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CEACAM5/yasara/00017567-Q137P.sce", download = "00017567-Q137P.sce", icon("download"), "Mutant CEACAM5 (Q137P) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CTNNA1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CTNNA1/yasara/00017566-R194K.sce", download = "00017566-R194K.sce", icon("download"), "Mutant CTNNA1 (R194K) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "CTNNB1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/CTNNB1/yasara/00017564-G34E.sce", download = "00017564-G34E.sce", icon("download"), "Mutant CTNNB1 (G34E) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "DCC"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/DCC/yasara/00017562-Q486L.sce", download = "00017562-Q486L.sce", icon("download"), "Mutant DCC (Q486L) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "EPCAM"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/EPCAM/yasara/00017570-I277M.sce", download = "00017570-I277M.sce", icon("download"), "Mutant EPCAM (I277M) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "FOS"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/FOS/yasara/00016603-V77G.sce", download = "00016603-V77G.sce", icon("download"), "Mutant FOS (V77G) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "KIT"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/KIT/yasara/00017574-D816V_PROB.sce", download = "00017574-D816V_PROB.sce", icon("download"), "Mutant KIT (D816V) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "KITLG"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/KITLG/yasara/00017579-D149Y.sce", download = "00017579-D149Y.sce", icon("download"), "Mutant KITLG (D149Y) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "KRAS"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/KRAS/yasara/00017581-G60R.sce", download = "00017581-G60R.sce", icon("download"), "Mutant KRAS (G60R) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "KRT20"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/KRT20/yasara/00016599-H208Y.sce", download = "00016599-H208Y.sce", icon("download"), "Mutant KRT20 (H208Y) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MALT1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MALT1/yasara/00017583-I60G.sce", download = "00017583-I60G.sce", icon("download"), "Mutant MALT1 (I60G) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MET"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MET/yasara/00017593-D1265Y.sce", download = "00017593-D1265Y.sce", icon("download"), "Mutant MET (D1265Y) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MGMT"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MGMT/yasara/00017591-R159Q.sce", download = "00017591-R159Q.sce", icon("download"), "Mutant MGMT (R159Q) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MMP2"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MMP2/yasara/00017589-C288F.sce", download = "00017589-C288F.sce", icon("download"), "Mutant MMP2 (C288F) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MSH2"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MSH2/yasara/00017638-G751R.sce", download = "00017638-G751R.sce", icon("download"), "Mutant MSH2 (G751R) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MTHFR"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MTHFR/yasara/00017600-R377C.sce", download = "00017600-R377C.sce", icon("download"), "Mutant MTHFR (R377C) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MUC1"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MUC1/yasara/00017606-V126M.sce", download = "00017606-V126M.sce", icon("download"), "Mutant MUC1 (V126M) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "MYC"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/MYC/yasara/00017604-P74A.sce", download = "00017604-P74A.sce", icon("download"), "Mutant MYC (P74A) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "PCNA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/PCNA/yasara/00017602-S39R.sce", download = "00017602-S39R.sce", icon("download"), "Mutant PCNA (S39R) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "PIK3CA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/PIK3CA/yasara/00016597-R555K.sce", download = "00016597-R555K.sce", icon("download"), "Mutant PIK3CA (R555K) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "PTEN"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/PTEN/yasara/00017636-D252G.sce", download = "00017636-D252G.sce", icon("download"), "Mutant PTEN (D252G) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "PTGS2"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/PTGS2/yasara/00017608-E488G.sce", download = "00017608-E488G.sce", icon("download"), "Mutant PTGS2 (E488G) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "RUNX3"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/RUNX3/yasara/00017613-G335R.sce", download = "00017613-G335R.sce", icon("download"), "Mutant RUNX3 (G335R) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "SDHA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/SDHA/yasara/00017611-R600Q.sce", download = "00017611-R600Q.sce", icon("download"), "Mutant SDHA (R600Q) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "SDHB"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/SDHB/yasara/00017609-C101Y.sce", download = "00017609-C101Y.sce", icon("download"), "Mutant SDHB (C101Y) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "SDHD"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/SDHD/yasara/00017628-H102L.sce", download = "00017628-H102L.sce", icon("download"), "Mutant SDHD (H102L) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "SMAD4"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/SMAD4/yasara/00017607-R361C.sce", download = "00017607-R361C.sce", icon("download"), "Mutant SMAD4 (R361C) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "STK11"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/STK11/yasara/00017625-D194N.sce", download = "00017625-D194N.sce", icon("download"), "Mutant STK11 (D194N) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "TNF"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/TNF/yasara/00017623-I194N.sce", download = "00017623-I194N.sce", icon("download"), "Mutant TNF (I194N) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "TP53"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/TP53/yasara/00017641-R282W.sce", download = "00017641-R282W.sce", icon("download"), "Mutant TP53 (R282W) .sce File")
                    
        )
      }
      else if(input$radio_func_yasara == "VEGFA"){
        tags$button(style="padding: 15px 25px;
  font-size: 15px;
  text-align: center;
  cursor: pointer;
  outline: none;
  color: #fff;
  background-color: navy;
  border: none;
  border-radius: 2px;
  box-shadow: 0 6px #999;",
                    
                    
                    a(href ="Functional/Free energy calculation/VEGFA/yasara/00017587-R288W.sce", download = "00017587-R288W.sce", icon("download"), "Mutant VEGFA (R288W) .sce File")
                    
        )
      }
      
      
      
    })
    
    
   
    
}

shinyApp(ui = ui, server = server)