function metaCycModel=getModelFromMetaCyc(metacycPath,keepTransportRxns,keepUnbalanced,keepUndetermined)
% getModelFromMetaCyc
%   Retrieves information stored in MetaCyc flat files and generates a super model
%
%
%   metacycPath         if metaCycRxns.mat is not in the RAVEN\external\metacyc
%                       directory, this function will attempt to read data
%                       from a local dump of the metacyc database. metacycPath
%                       is the path to the root of this database
%   keepTransportRxns   include transportation reactions, which often have identical
%                       reactants and products that turn to be all-zero columns in
%                       the S matrix (opt, default false)
%   keepUnbalanced      include reactions cannot be unbalanced reactions, usually
%                       because they are polymeric reactions or because of a
%                       specific difficulty in balancing class structures
%                       (opt, default false)
%   keepUndetermined    include reactions that have substrates lack chemical
%                       structures or with non-numerical coefficients (e.g. n+1)
%                       (opt, default false)
%
%   metaCycModel        a model structure generated from MetaCyc database
%                       including all reactions and metabolites in MetaCyc
%
%   Usage: getModelFromMetaCyc(metacycPath)
%
%   Hao Wang, 2017-06-17
%

if nargin<1
    metacycPath='';
end
if nargin<2
    keepTransportRxns=false;
end
if nargin<3
    keepUnbalanced=false;
end
if nargin<4
    keepUndetermined=false;
end

%First get all reactions
metaCycModel=getRxnsFromMetaCyc(metacycPath,keepTransportRxns,keepUnbalanced,keepUndetermined);
fprintf('MetaCyc reactions loaded\n');

%Get reaction and enzyme association
metaCycEnzymes=getEnzymesFromMetaCyc(metacycPath);
fprintf('MetaCyc enzymes loaded\n');

%Replace rxnNames with those from metaCycEnzymes
[a, b]=ismember(metaCycModel.rxns,metaCycEnzymes.rxns);
a=find(a);
b=b(a);
metaCycModel.rxnNames(a)=metaCycEnzymes.rxnNames(b);

%Create the rxnGeneMat for the reactions, by geting all enzymes and
%corresponding subunits
rxnNum=numel(metaCycModel.rxns);
metaCycModel.genes=metaCycEnzymes.enzymes;
metaCycModel.rxnGeneMat=sparse(rxnNum,numel(metaCycEnzymes.enzymes)); % row: rxn, column: enzyme
metaCycModel.grRules=cell(rxnNum,1);

%Loop through all reactions to generate rxnGeneMat matrix and grRules
%This step also cross-link reactions to their catalyzing enzymes
for i=1:rxnNum

	metaCycModel.grRules{i}='';
	%Find out if this is an enzymatic reaction
	[a b]=ismember(metaCycModel.rxns(i),metaCycEnzymes.rxns);
	if a
		I=[];   %Find out all catalyzing enzymes, which are treated as isoenzymes
		I=find(metaCycEnzymes.rxnEnzymeMat(b,:));
		if ~isempty(I)
			
			grRule='';
			for j=1:numel(I)
								
				subgrRule=''; %Find out if enzyme complex
				[c d]=ismember(metaCycEnzymes.enzymes(I(j)),metaCycEnzymes.cplxs);
				if c
						for k=1:numel(metaCycEnzymes.cplxComp{d}.subunit)
							if strcmp(subgrRule,'')
								subgrRule=strcat('(',metaCycEnzymes.cplxComp{d}.subunit{k});
							else
								subgrRule=strcat(subgrRule,{' and '},metaCycEnzymes.cplxComp{d}.subunit{k});
							end
							
							[x geneIndex]=ismember(metaCycEnzymes.cplxComp{d}.subunit{k},metaCycModel.genes);
							%if x
								metaCycModel.rxnGeneMat(i,geneIndex)=1;
							%else
							%	disp(metaCycEnzymes.cplxComp{d}.subunit{k});
							%end
						end
						subgrRule=strcat(subgrRule,')');
				else
						subgrRule=metaCycEnzymes.enzymes(I(j));
						metaCycModel.rxnGeneMat(i,I(j))=1;
				end
				
				%Generating grRules
				if ~strcmp(subgrRule,'')
					if ~strcmp(grRule,'')
						grRule=strcat(grRule,{' or '},subgrRule);
					else
						grRule=subgrRule;
					end
				end
				
			end
			metaCycModel.grRules{i}=grRule;
						
		end

	end
end

%Then get all metabolites
metaCycMets=getMetsFromMetaCyc(metacycPath);
fprintf('MetaCyc compounds loaded\n');

%Add information about all metabolites to the model
[a, b]=ismember(metaCycModel.mets,metaCycMets.mets);
a=find(a);
b=b(a);

if ~isfield(metaCycModel,'metNames')
   metaCycModel.metNames=cell(numel(metaCycModel.mets),1);
   metaCycModel.metNames(:)={''};
end
metaCycModel.metNames(a)=metaCycMets.metNames(b);

if ~isfield(metaCycModel,'metFormulas')
   metaCycModel.metFormulas=cell(numel(metaCycModel.mets),1);
   metaCycModel.metFormulas(:)={''};
end
metaCycModel.metFormulas(a)=metaCycMets.metFormulas(b);

if ~isfield(metaCycModel,'metCharge')
   metaCycModel.metCharge=zeros(numel(metaCycModel.mets),1);
end
metaCycModel.metCharge(a)=metaCycMets.metCharge(b);

if ~isfield(metaCycModel,'inchis')
   metaCycModel.inchis=cell(numel(metaCycModel.mets),1);
   metaCycModel.inchis(:)={''};
end
metaCycModel.inchis(a)=metaCycMets.inchis(b);

if ~isfield(metaCycModel,'metMiriams')
   metaCycModel.metMiriams=cell(numel(metaCycModel.mets),1);
end
metaCycModel.metMiriams(a)=metaCycMets.metMiriams(b);

if ~isfield(metaCycModel,'keggid')
   metaCycModel.keggid=cell(numel(metaCycModel.mets),1);
end
metaCycModel.keggid(a)=metaCycMets.keggid(b);

%Put all metabolites in one compartment called 's' (for system). This is
%done just to be more compatible with the rest of the code
metaCycModel.comps={'s'};
metaCycModel.compNames={'System'};
metaCycModel.compOutside={''};
metaCycModel.metComps=ones(numel(metaCycModel.mets),1);

%If reactions with undefined stoichiometry are kept, then the corresponding
%metabolites will have ids such as "(n+1) C000001" and their names will be
%empty. These ids are not valid SBML identifiers and are therefore replaced
%with "undefined1, undefined2...". The former ids are stored as the new
%names
%I=find(cellfun(@any,strfind(metaCycModel.mets,'n')) | cellfun(@any,strfind(metaCycModel.mets,'m')));
%I=find(cellfun(@any,strfind(metaCycModel.mets,'n')) | cellfun(@any,strfind(metaCycModel.mets,'m')) | cellfun(@any,strfind(metaCycModel.mets,'n+1')) | cellfun(@any,strfind(metaCycModel.mets,'m+1')));
%metaCycModel.metNames(I)=metaCycModel.mets(I);
%repNums=1:numel(I);
%repIDs=strcat('undefined_',cellfun(@num2str,num2cell(repNums(:)),'UniformOutput',false));
%metaCycModel.mets(I)=repIDs;

%It could also be that the metabolite and reaction names are empty
%for some reasons. In that case, use the ID instead
I=cellfun(@isempty,metaCycModel.metNames);
metaCycModel.metNames(I)=metaCycModel.mets(I);
I=cellfun(@isempty,metaCycModel.rxnNames);
metaCycModel.rxnNames(I)=metaCycModel.rxns(I);

end