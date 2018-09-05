%This exercise is to reconstruct draft GEMs by using both KEGG and MetaCyc
%pathway databases. A combined model with comprehensive coverage of
%metabolic pathways will be eventually generated from different de novo
%reconstruction approaches. The input is a FASTA format file with
%whole-proteome sequences. The combined model will be subsequently used
%for refinement of existing high-quality model and generation of a new
%version of GEM, by utilizing the manual curations results. This tutorial
%is a showcase of the new features released in REVAN 2.0 through 
%demonstrating the utilization of the newly developed functions on GEM
%reconstruction and curation for Streptomyces coelicolor strain A3(2).
%Uers may apply the introduced approaches in their own work for other
%organisms.
%
% Hao Wang, 2018-08-15
%

%Noted that the detail instructions of RAVEN functions can be found both
%by the Matlab command 'help functionName' and from opening the HTML file
%functionName.html under the 'doc' folder. The available parameters,
%input and output files, the choic of parameters and their default values
%are specified in the function descriptions. Any futher questions can be
%brought up through Gitter chatting room or posted as an issues on GitHub
%from: https://github.com/SysBioChalmers/RAVEN/issues/new

%Before reconstruction, a FASTA file with protein sequences of the target
%organism need to be prepared. In this tutorial, all protein sequences
%of S. coelicolor A3(2) were downloaded from the NCBI genome database and
%provided in current folder (Sco_all_protein.faa). The description lines
%(starting with ">" character) of this FASTA file were modified by
%keeping only the locus tag (e.g. SCO6005), which is commonly used as
%the gene identifiers in GEMs.

%This command generates a draft GEM using the MetaCyc database. The first
%parameter is organismID and need to be speicified by user. The other
%parameters are set to exclude unbalanced and undetermined reactions,
%but keep transport rections. The two parameters for homology search are
%set to default values that have been optimized to capture the protein
%hits with both comprehensive coverage and the least false positives.
ScoMetaCycDraftModel=getMetaCycModelForOrganism('ScoMetaCyc','Sco_all_protein.faa',1);

%Two draft models are generated by using the KEGG database. The first one
%is reconstructed from the genome information annotated by KEGG. Since
%KEGG provides genome annotations for over 5000 organisms, their GEMs
%can thus be generated without providing protein sequnces and with only
%the organismID. The organisms with KEGG annotation and their associated
%ids (e.g. 'sco' for Streptomyces coelicolor) can be found from here:
%https://www.kegg.jp/kegg/catalog/org_list.html
%The following command generates a draft model using KEGG annotation and
%exclude incomplete reactions and reactions with undefined stoichiometry.
ScoKEGGAnnotation=getKEGGModelForOrganism('sco','','','',0,0);

%The second KEGG-based draft model is generated based on sequence homology
%to KEGG Ortholog sequence clusters with excluding imcomplete reactions and
%the ones with undefined stoichiometry. Type "help getKEGGModelForOrganism"
%to see the detail instructions for the choice of different parameters.
%The default values for homology search are used because they have been
% optimized for the best performance.
ScoKEGGHomology=getKEGGModelForOrganism('ScoKEGGHMMs','Sco_all_protein.faa','prok90_kegg82','',0,0);

%Each of the three de novo reconstruction process takes about 10 minutes
%on a regular computer. Given that KEGG and MetaCyc databases are formulated
%in different ways; KEGG relys on sequence-based annotaiton, while MetaCyc
%collects only experimentally verified pathways. Therefore, integration of
%MetaCyc- and KEGG-based draft models could have a better coverage of the
%metbolism for the target organism.

%At first, the two KEGG-based models can be directly merged
ScoKEGGDraftModel=mergeModels({ScoKEGGAnnotation ScoKEGGHomology});

numel(ScoKEGGHomology.rxns)+numel(ScoKEGGAnnotation.rxns)
numel(ScoKEGGDraftModel.rxns)
%By checking the reaction number, it can be found that reaction number
%in the merged model equals adding up the reaction numbers in homology
%and annotation KEGG draft models. And the there are duplicated reactions
%in this merged model.

%To merge the duplicated reactions into one and combine multiple
%iso-enzymes into a single grRules. The expandModel and contractModel
%functions are applied, see their instructions for detail.
ScoKEGGDraftModel=expandModel(ScoKEGGDraftModel);
ScoKEGGDraftModel=contractModel(ScoKEGGDraftModel);

%In the end, KEGG- and MetaCyc-based draft models can be further combined
%into an integrated GEM with 2605 reactions, 3005 metabolites and 2175 genes.
%This step is achieved by the function combineMetaCycKEGGModels, which
%fristly converts metabolite and reactions identifiers in the KEGG model
%into corresponding MetaCyc ids, and then detect duplications and keep
%only unique reactions and metabolites that are mostly in MetaCyc
%identifiers.
ScoCombinedDraftModel=combineMetaCycKEGGModels(ScoMetaCycDraftModel, ScoKEGGDraftModel)

%With this combined model, we are going to refine an existing high
%quality GEM iMK1208 by incorporating the new pathways/reactions that
%are found in the combined model but absent from the previous iMK1208
%model. A total of 398 reactions in the combined draft model were
%determined as new pathways based on manual curation results, which
%have been organized into the Excel file SupportingTables.xlsx.
%Now read in these manually selected reactions and their subSystems
%from the sheet TableS3 into an array structure selectedNewRxns
[~, textData]=xlsread('SupportingTables.xlsx','TableS3');
selectedNewRxns.rxns=textData(2:end,1);
selectedNewRxns.subSystems=textData(2:end,3);
selectedNewRxns

%Next step is to generate a sub-model that includes only these new
%reactions. This is implemented by subtracting the other reactions
%from the combined model using function removeReactions
rxnsToRemove=setdiff(ScoCombinedDraftModel.rxns,selectedNewRxns.rxns);
newRxnSubModel=removeReactions(ScoCombinedDraftModel,rxnsToRemove,1,1);

%Now this newRxnSubModel contains only these new reactions. It needs
%to be modified before merging with the iMK1208. Since these reactions
%are metabolic ones and can all be assigned to the cytoplasm compartment
newRxnSubModel.comps{1}='c';
newRxnSubModel.compNames{1}='Cytoplasm';

%Here are some amendaments to the genes and rxnGeneMat fields
newRxnSubModel.genes={};
for k=1:numel(newRxnSubModel.rxns)
   newRxnSubModel.genes=[newRxnSubModel.genes;transpose(strsplit(newRxnSubModel.grRules{k},' or '))];
end
newRxnSubModel.genes=unique(newRxnSubModel.genes);
[~, newRxnSubModel.rxnGeneMat, ~]=standardizeGrRules(newRxnSubModel);

%Since the newRxnSubModel is ready for incorportation, the iMK1208 model
%can be read in for integration
load('iMK1208.mat');

%It should be noted that the incompatible nomenclatures (especially the
%metabolite identifiers) used in different GEMs and databases led to a 
%serious problme in model comparison, curation and integration. In order
%to properly integrate the iMK1208 with this sub-modle. Both models
%should be unified into the same name space. Since the metabolites in the
%sub-model use MetaCyc identifiers that are different from the ones used
%in iMK1208 model. We implemented database mining and intensive manual
%curation in associating these metabolite identifiers and organized the
%results into sheet TableS2 in the SupportingTables.xlsx.
[~, textData]=xlsread('SupportingTables.xlsx','TableS2');
metaCycMetsIniMK=textData(2:end,4);

%The following section is to replace the metabolites in the sub-model
%with the identifiers and names used in iMK1208 according to the mapping
%information in TableS2.
[a, b]=ismember(newRxnSubModel.mets,metaCycMetsIniMK);
I=find(a);
newRxnSubModel.mets=strcat(newRxnSubModel.mets,'_c');
newRxnSubModel.metNames(I)=iMK1208.metNames(b(I));
newRxnSubModel.mets(I)=iMK1208.mets(b(I));


%We have resolved the necessary issues for merging of iMK1208 with the
%new reactions determed from de novo reconstructions. They can be
%direclty merged into the Sco4 model, which refers to the 4th major
%update of the GEM for Streptomyces coelicolor.
Sco4=mergeModels({iMK1208 newRxnSubModel});

sol=solveLP(Sco4)
%Let's check out if the newly generated Sco4 could grow or not. And it
%turn out this new model is functional.


%Here in this tutorial we present a complete reconstruciton process using
%different approaches introduced by RAVEN 2.0. Furthermore, the combined draft
%model obtained from de novo reconstructions is used in model refinments
%and curation for making an upgraded version of GEM.

