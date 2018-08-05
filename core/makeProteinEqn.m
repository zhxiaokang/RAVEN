function eqn = makeProteinEqn(proteinName,inputSeq)
% makeProteinEqn
%   Parses protein forming equation for proteinName using inputSeq
%   
%   proteinName     String containing protein name
%   inputSeq        String containing amino acid sequence of protein
%
%   Example: If multiple proteins are being added
%
%   protStruct = {}; %create empty protein structure
%   protStruct.name = {'MFalpha2','Myrosinase','P28'}; %add protein names
%   protStruct.seq = {'MKFISTFLTFILAAVSVTASSDEDIAQVPAEAIIGYLDFGGDHDIAFLPFSNATASGLLFINTTIAEAAEKEQNTTLAKREAVADAWHWLNLRPGQPMYKREANADAWHWLQLKPGQPMY', ...
%    'MKHLGLILAFLLALATCKADEEITCEENLPFKCSQPDRLNSSSFEKDFIFGVASSAYQACCLGRGLNVWDGFTHRYPNKSGPDHGNGDTTCDSFSYWQKDIDVLDELNATGYRFSIAWSRIIPRGKRSRGVNKDGINYYHGLIDGLIDKGITPFVTLFHWDLPQVLQDEYEGFLDPQIIHDFKHYANLCFQEFGHKVKNWLTINQLYTVPTRGYGAGSDAPGRCSPMVDPTCYAGNSSTEPYIVAHNQLLAHATVVDLYRKNYSIGPVMITRWFLPYNDTDPDSIAATERMKEFFLGWFMGPLTNGTYPQIMIDTVGERLPSFSPEESNLVKGSYDYLGLNYYVTQYAQPSPNPVHWANHTAMMDAGAKLTFRGNSDETKNSYYYPKGIYYVMDYFKTKYYNPLIYVTENGISTPGNETRDESMLHYKRIEYLCSHLCFLSKVIKEKHVNVKGYFAWSLGDNYEFDKGFTVRFGLSYIDWNNVTDRDLKLSGKWYQKFISPAIKNPLKKDFLRSSLTFEKNKKFEDA', ...
%    'LSTAADMQGVVTDGMASGLDKDYLKPDD'}; %add amino acid sequences
%
%   for x = 1:length(protStruct.name) %create an equation for each protein using makeProteinEqn()
%      protStruct.eqn{x} = makeProteinEqn(protStruct.name{x},protStruct.seq{x});
%   end
%
%   Usage: eqn = makeProteinEqn(proteinName,inputSeq)
%
% NOTE: This function was written for the metabolite nomenclature of the
% Chalmers Sysbio Saccharomyces cereviae GEM available here:
% https://github.com/SysBioChalmers/yeast-GEM. The amino acid structure
% (AAstruct) should be modified accordingly if a different amino acid 
% nomenclature is being used.
%
% Francisco Zorrilla, 2018-08-02 

AAstruct = {'Ala-tRNA(Ala)[c]' 'A'
'Arg-tRNA(Arg)[c]' 'R'
'Asn-tRNA(Asn)[c]' 'N'
'Asp-tRNA(Asp)[c]' 'D'
'Cys-tRNA(Cys)[c]' 'C'
'Gln-tRNA(Gln)[c]' 'Q'
'Glu-tRNA(Glu)[c]' 'E'
'Gly-tRNA(Gly)[c]' 'G'
'His-tRNA(His)[c]' 'H'
'Ile-tRNA(Ile)[c]' 'I'
'Leu-tRNA(Leu)[c]' 'L'
'Lys-tRNA(Lys)[c]' 'K'
'Met-tRNA(Met)[c]' 'M'
'Phe-tRNA(Phe)[c]' 'F'
'Pro-tRNA(Pro)[c]' 'P'
'Ser-tRNA(Ser)[c]' 'S'
'Thr-tRNA(Thr)[c]' 'T'
'Trp-tRNA(Trp)[c]' 'W'
'Tyr-tRNA(Tyr)[c]' 'Y'
'Val-tRNA(Val)[c]' 'V'
};

for q = 1:length(AAstruct) %create third column with all zeros
    AAstruct{q,3}=0;
end
    
ATPcost = num2str((length(inputSeq)*4)-1);

for c = 1:length(inputSeq) %cycle through each AA in input sequence
    for d = 1:length(AAstruct) %cycle through each possible AA
       if inputSeq(c) == AAstruct{d,2} %check if match
           AAstruct{d,3}= AAstruct{d,3}+1; %add count if match
       end
    end
end

%for e = 1:length(AAstruct)   ATTEMPT FILTERING OUT AAs WITH 0 COUNTS,
%                             OTHERWISE MANUALLY PULL THEM OUT LATER
%   if AAstruct{e,3} == 0
%       AAstruct{e,1} = {}
%       AAstruct{e,2} = {}
%       AAstruct{e,3} = {}
%   end
%end

for w = 1:length(AAstruct) %convert numbers to strings
    AAstruct{w,3} = num2str(AAstruct{w,3});
end

rxtntSide = ''; %initialize vectors to hold eqn
prodSide = ''; 
eqn = ''; 

for r = 1:length(AAstruct) % generate rxtns side of eqation
    rxtntSide = [rxtntSide ' + ' AAstruct{r,3} ' ' AAstruct{r,1}];
    prodSide = [prodSide ' + ' AAstruct{r,3} ' ' AAstruct{r,1}(5:length(AAstruct{r,1}))];
end

eqn = [rxtntSide(4:length(rxtntSide)) ' + ' ATPcost ' ATP[c]' ' + ' ATPcost ' H2O[c] => ' prodSide(4:length(prodSide)) ' + ' ATPcost ' ADP[c]' ' + ' ATPcost ' phosphate[c]' ' + ' ATPcost ' H+[c]' ' + ' proteinName '[c]' ]; %add ATP

end