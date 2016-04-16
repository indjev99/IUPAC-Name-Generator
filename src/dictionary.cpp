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
    curr_dict.help+="\nUse the Right mouse button to cancel the current action and remove atoms.\nHold down O, F, C, B or I in order to place aн Oxygen, Fluorine, Chlorine,\nBromine or Iodine atom respectively.";
    curr_dict.help+="\nPress Backspace to undo.\nPress Shift to toggle snapping on and off.\nPress R to reset.\nPress L to change the language.";
    curr_dict.help+="\nPress N to enter a name of a compound::\nPress H for help.\nPress Escape to end the program.";
    curr_dict.PACTC="Press any key to continue.";
    dictionaries.push_back(curr_dict);

    curr_dict.CNP= {"алка","мета","ета","пропа","бута"};
    curr_dict.CBI= {curr_dict.error,"ан","ен","ин"};
    curr_dict.NP= {curr_dict.error,"моно","ди","три","тетра","пента","хекса","хепта","окта","нона",
                   "дека","ундека","додека","тридека","тетрадека","пентадека","хексадека","хептадека","октадека","нонадека"
                  };
    curr_dict.SS= {curr_dict.error,"ил","илиден","илиден"};
    curr_dict.CP="цикло";

    curr_dict.HP["F"]="флуоро";
    curr_dict.HP["Cl"]="хлоро";
    curr_dict.HP["Br"]="бромо";
    curr_dict.HP["I"]="йодо";

    curr_dict.FGS["alcohol"]="ол";
    curr_dict.FGP["alcohol"]="хидрокси";
    curr_dict.FGS["ketone"]="он";
    curr_dict.FGP["ketone"]="оксо";
    curr_dict.FGS["aldehyde"]="ал";
    curr_dict.FGP["aldehyde"]="оксо";
    curr_dict.FGS["carboxylic acid"]="ова киселина";
    curr_dict.FGP["carboxylic acid"]="карбокси";

    curr_dict.benzene="бензен";
    curr_dict.phen="фен";

    curr_dict.NC="Не са свързани";
    curr_dict.help="Използвайте средния бутон на мишката, за да започнете или продължите вериги.\nИзползвайте левия бутон на мишката, за да завършите вериги\nили да местите атоми.";
    curr_dict.help+="\nИзплозвайте десния бутон на мишката, за да откажете текущото действие\nили да махнете атом.\nЗадръжте O, F, C, B или I, за да сложите кислороден, флуорен, хлорен,\nбромен или йоден атом съответно.";
    curr_dict.help+="\nНатиснете Backspace, за да върнете последното действие.\nНатиснете Shift, за да включите или изключите\nавтоматичното наместване на атомите.\nНатиснете R, за да рестартирате.";
    curr_dict.help+="\nНатиснете L, за да смените езика.\nНатиснете H за инструкции.\nНатиснете N, за да въведете името на съединение.";
    curr_dict.help+="\nНатиснете Escape, за да изключите програмата.";
    curr_dict.PACTC="Натиснете който и да е клавиш, за да продължите.";
    dictionaries.push_back(curr_dict);

    curr_dict=dictionaries[0];
}
