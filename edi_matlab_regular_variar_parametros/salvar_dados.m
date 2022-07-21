function []=salvar_dados(cell_dados)
%% Abrir planilha para salvar os dados
planilha = 'analise_decolagem.xlsx';
aba_planilha = 'main';
%% Numero de analises feitas variando a deflexao do profundor
n_de = size(cell_dados,1);
%% Escrever os valores de ac_eh e massa, separados por deflexao
col_num = 1;
for i=1:n_de
    %% Numero de analises feitas variando a distancia de acionamento
    n_ac_eh = size(cell_dados{i,1},1);
    %% Transformar índice em endereço A1 do Excel
    col_address_i=num2xlcol(col_num);
    col_address_f=num2xlcol(col_num+2);
    linha_inicial = int2str(2);
    linha_final = int2str(n_ac_eh+1);
    xlRange = strcat(col_address_i,linha_inicial,':',col_address_f,linha_final);
    %% Escrever ac_eh X MTOW
    matriz_output = [cell_dados{i,1}(:,1) cell_dados{i,1}(:,2) cell_dados{i,1}(:,3)];
    xlswrite(planilha,matriz_output,aba_planilha,xlRange);
    %% Escrever deflexao
    de_output = cell_dados{i,2};
    xlRange = strcat(col_address_i,linha_inicial-1);
    xlswrite(planilha,de_output,aba_planilha,xlRange);
    col_num = col_num + 4;
end
end