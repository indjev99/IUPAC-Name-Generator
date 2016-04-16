#include "dictionary.h"
string dictionary::getCNP(int CN)
{
    if (CN>=CNP.size())
    {
        if (CN>=NP.size()) CN=0;
        else return NP[CN];
    }
    return CNP[CN];
}
string dictionary::getNP(int ST)
{
    if (ST>=NP.size()) ST=0;
    return NP[ST];
}
void setDictionaries()
{
    curr_dict.CNP= {"alka","metha","etha","propa","buta"};
    curr_dict.CBI= {curr_dict.error,"ane","ene","yne"};
    curr_dict.NP= {curr_dict.error,"mono","di","tri","tetra","penta","hexa","hepta","octa","nona",
                   "deca","undeca","dodeca","trideca","tetradeca","pentadeca","hexadeca","heptadeca","octadeca","nonadeca",
                   "icosa","henicosa","docosa","tricosa","tetracosa","pentacosa","hexacosa","heptacosa","octacosa","nonacosa",
                   "triaconta","hentriaconta","hentriaconta","tritriaconta"
                  };
    curr_dict.SS= {curr_dict.error,"yl","ylidene","ylidyne"};
    curr_dict.CP="cyclo";

    curr_dict.HP["F"]="fluoro";
    curr_dict.HP["Cl"]="chloro";
    curr_dict.HP["Br"]="bromo";
    curr_dict.HP["I"]="iodo";

    curr_dict.FGS["alcohol"]="ol";
    curr_dict.FGP["alcohol"]="hydroxy";
    curr_dict.FGS["ketone"]="one";
    curr_dict.FGP["ketone"]="oxo";
    curr_dict.FGS["aldehyde"]="al";
    curr_dict.FGP["aldehyde"]="oxo";
    curr_dict.FGS["carboxylic acid"]="oic acid";
    curr_dict.FGP["carboxylic acid"]="carboxy";

    curr_dict.benzene="benzene";
    curr_dict.phen="phen";

    curr_dict.NC="Not Connected";
    curr_dict.help="Use the Middle Mouse Button to start or continue chains.\nUse the Left Mouse Button to end chains and to move atoms.";
    curr_dict.help+="\nUse the Right mouse button to cancel the current action and remove atoms.\nHold down O, F, C, B or I in order to place a� Oxygen, Fluorine, Chlorine,\nBromine or Iodine atom respectively.";
    curr_dict.help+="\nPress Backspace to undo.\nPress Shift to toggle snapping on and off.\nPress R to reset.\nPress L to change the language.";
    curr_dict.help+="\nPress N to enter a name of a compound::\nPress H for help.\nPress Escape to end the program.";
    curr_dict.PACTC="Press any key to continue.";
    dictionaries.push_back(curr_dict);

    curr_dict.CNP= {"����","����","���","�����","����"};
    curr_dict.CBI= {curr_dict.error,"��","��","��"};
    curr_dict.NP= {curr_dict.error,"����","��","���","�����","�����","�����","�����","����","����",
                   "����","������","������","�������","���������","���������","���������","���������","��������","��������"
                  };
    curr_dict.SS= {curr_dict.error,"��","������","������"};
    curr_dict.CP="�����";

    curr_dict.HP["F"]="������";
    curr_dict.HP["Cl"]="�����";
    curr_dict.HP["Br"]="�����";
    curr_dict.HP["I"]="����";

    curr_dict.FGS["alcohol"]="��";
    curr_dict.FGP["alcohol"]="��������";
    curr_dict.FGS["ketone"]="��";
    curr_dict.FGP["ketone"]="����";
    curr_dict.FGS["aldehyde"]="��";
    curr_dict.FGP["aldehyde"]="����";
    curr_dict.FGS["carboxylic acid"]="��� ��������";
    curr_dict.FGP["carboxylic acid"]="��������";

    curr_dict.benzene="������";
    curr_dict.phen="���";

    curr_dict.NC="�� �� ��������";
    curr_dict.help="����������� ������� ����� �� �������, �� �� ��������� ��� ���������� ������.\n����������� ����� ����� �� �������, �� �� ��������� ������\n��� �� ������� �����.";
    curr_dict.help+="\n����������� ������ ����� �� �������, �� �� �������� �������� ��������\n��� �� ������� ����.\n�������� O, F, C, B ��� I, �� �� ������� ����������, �������, ������,\n������ ��� ����� ���� ���������.";
    curr_dict.help+="\n��������� Backspace, �� �� ������� ���������� ��������.\n��������� Shift, �� �� �������� ��� ���������\n������������� ���������� �� �������.\n��������� R, �� �� ������������.";
    curr_dict.help+="\n��������� L, �� �� ������� �����.\n��������� H �� ����������.\n��������� N, �� �� �������� ����� �� ����������.";
    curr_dict.help+="\n��������� Escape, �� �� ��������� ����������.";
    curr_dict.PACTC="��������� ����� � �� � ������, �� �� ����������.";
    dictionaries.push_back(curr_dict);

    curr_dict=dictionaries[0];
}
