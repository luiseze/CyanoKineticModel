function model = MJanasch_MDF2_LoadData(fileName)

%% Read in the excel file
data_raw = importdata(fileName);

%% Create Model Structure
% Reactions
model.reactions.rxns        = data_raw.textdata.Reactions(2:end,2);
model.reactions.equations   = data_raw.textdata.Reactions(2:end,4);
model.reactions.lb          = [];
model.reactions.ub          = [];
model.reactions.dG0         = data_raw.data.Reactions;

%Rest potentially added later

% Metabolites
model.metabolites.mets      = data_raw.textdata.Metabolites(2:end,2);

model.metabolites.lconc     = data_raw.data.Metabolites(1:end,1);
model.metabolites.uconc     = data_raw.data.Metabolites(1:end,2);

%% Create S-matrix
% Initialize the stoichiometric matrix
reaction_string = data_raw.textdata.Reactions(2:end,4);
mets            = data_raw.textdata.Metabolites(2:end,2);
rxns            = model.reactions.rxns;

S = zeros(numel(mets),numel(reaction_string));
reaction_string=strrep(reaction_string,'+', '¤');
    
% Loop through the reaction_string and add the info to the S matrix
for i = 1:numel(reaction_string)
    % Start by finding the position of the (=> or <=>)
    arrowIndex = strfind(reaction_string{i},'<=>');
    
    % Split reactions into reactants and products
    substrates  = regexp(reaction_string{i}(1:arrowIndex-1),'¤','split');
    products    = regexp(reaction_string{i}(arrowIndex+3:end),'¤','split');
    
    % If the splitting character is at the end (if exchange rxns), then an
    % empty string will exist together with the real ones. Remove it
    substrates(cellfun(@isempty,substrates))=[];
    products(cellfun(@isempty,products))=[];

    % A vector where an element is -1 is the corresponding metabolite is a
    % reactant and 1 if it's a product
    multiplyWith=[ones(numel(substrates),1)*-1; ones(numel(products),1)];

    metabolites=[substrates products];

    % Now loop through the metabolites and see if the metabolite has a coefficient 
    % (it will look as 'number name')
        for j=1:numel(metabolites)
            space=strfind(metabolites{j},' ');
        
            if isempty(space)
                %No coefficient
                coeff=1;
                name=metabolites{j};
            else
                coeff=str2double(metabolites{j}(1:space(1)));
            
                %If it was not a coefficiant
                if isnan(coeff)
                    coeff=1;
                    name=metabolites{j};
                else
                    name=metabolites{j}(space+1:end);
                end
            end
        
            %Find the name in the mets list
            %[a b]=ismember(name,mets);
            b=find(strcmp(name,mets),1);
        
            if any(b)
                S(b,i)=S(b,i)+coeff*multiplyWith(j);
            else
                if isempty(rxns)
                    dispEM(['Could not find metabolite ' name ' in metabolite list']);
                else
                    dispEM(['The metabolite "' name '" in reaction ' rxns{i} ' was not found in the metabolite list']);
                end
            end
        end
end
S=sparse(S);

% Add S-matrix to model structure
model.rawS = S;

end
