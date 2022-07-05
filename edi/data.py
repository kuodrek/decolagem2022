import openpyxl as xl

def get_data():
    workbook = xl.load_workbook(filename="dadoscompletosMONOPLANO.xlsx", data_only=True)
    vars_dados = workbook['vars']

    # Variáveis globais
    timestep = vars_dados.cell(row=3,column=4).value
    headwind = vars_dados.cell(row=4,column=4).value
    TOW = vars_dados.cell(row=5,column=4).value
    n_helice = vars_dados.cell(row=6,column=4).value
    altitude = vars_dados.cell(row=7,column=4).value
    mi = vars_dados.cell(row=8,column=4).value
    Sd = vars_dados.cell(row=9,column=4).value
    ac_eh = vars_dados.cell(row=10,column=4).value

    # Variáveis da asa
    Sw = vars_dados.cell(row=3,column=5).value
    iw = vars_dados.cell(row=3,column=6).value
    CLmax_w = vars_dados.cell(row=3,column=7).value
    CL_0_w = vars_dados.cell(row=3,column=8).value
    CL_alfa_w = vars_dados.cell(row=3,column=9).value
    CD_fus = vars_dados.cell(row=3,column=10).value
    Cm_ac_w = vars_dados.cell(row=3,column=11).value
    MAC_w = vars_dados.cell(row=3,column=12).value
    ## Dados para polyfit das curvas de arrasto parasita e induzido
    v1 = vars_dados.cell(row=5,column=5).value
    v2 = vars_dados.cell(row=6,column=5).value
    CDp_w1 = vars_dados.cell(row=5,column=6).value
    CDp_w2 = vars_dados.cell(row=6,column=6).value

    a_1 = vars_dados.cell(row=5,column=7).value
    a_2 = vars_dados.cell(row=6,column=7).value
    a_3 = vars_dados.cell(row=7,column=7).value
    a_4 = vars_dados.cell(row=8,column=7).value
    a_5 = vars_dados.cell(row=9,column=7).value
    a_6 = vars_dados.cell(row=10,column=7).value

    CDi_1 = vars_dados.cell(row=5,column=8).value
    CDi_2 = vars_dados.cell(row=6,column=8).value
    CDi_3 = vars_dados.cell(row=7,column=8).value
    CDi_4 = vars_dados.cell(row=8,column=8).value
    CDi_5 = vars_dados.cell(row=9,column=8).value
    CDi_6 = vars_dados.cell(row=10,column=8).value

get_data()