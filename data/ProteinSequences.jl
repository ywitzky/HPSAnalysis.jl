module ProteinSequences

NameToSeq=Dict("RS31"=>"MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWAKGERGRPRGDAKAPSNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYALKDDDERDDRNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRARSPEYDRYKGPAAYERRRSPDYGRRSSDYGRQRSPGYDRYRSRSPVPRGRP",
"RS31a"=>"MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWAKDFQGERGKPRDGKAVSNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYALREAGEREDRYAGSRRRRSPSPVYRRRPSPDYTRRRSPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG",
"RS40"=>"MKPVFCGNFEYDAREGDLERLFRKYGKVERVDMKAGFAFVYMEDERDAEDAIRALDRFEFGRKGRRLRVEWTKSERGGDKRSGGGSRRSSSSMRPSKTLFVINFDADNTRTRDLEKHFEPYGKIVNVRIRRNFAFIQYEAQEDATRALDASNNSKLMDKVISVEYAVKDDDARGNGHSPERRRDRSPERRRRSPSPYKRERGSPDYGRGASPVAAYRKERTSPDYGRRRSPSPYKKSRRGSPEYGRDRRGNDSPRRRERVASPTKYSRSPNNKRERMSPNHSPFKKESPRNGVGEVESPIERRERSRSSPENGQVESPGSIGRRDSDGGYDGAESPMQKSRSPRSPPADE",
"RS41"=>"MKPVFCGNFEYDARESDLERLFRKYGKVERVDMKAGFAFVYMEDERDAEDAIRALDRFEYGRTGRRLRVEWTKNDRGGAGRSGGSRRSSSGLRPSKTLFVINFDAQNTRTRDLERHFEPYGKIVNVRIRRNFAFIQYEAQEDATRALDATNSSKLMDKVISVEYAVKDDDSRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESSPSRSPPAEE",
"RS41_PHOS6"=>"MKPVFCGNFEYDARESDLERLFRKYGKVERVDMKAGFAFVYMEDERDAEDAIRALDRFEYGRTGRRLRVEWTKNDRGGAGRSGGSRRSSSGLRPSKTLFVINFDAQNTRTRDLERHFEPYGKIVNVRIRRNFAFIQYEAQEDATRALDATNSSKLMDKVISVEYAVKDDDSRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGA#PVAHKRER&#PDYGRGRRSPSPYKRARL#PDYKRDDRRRERVA#PENGAVRNR#PRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESSPSRSPPAEE",
"A1_LCD"=>"GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQSSNFGPMKGGNFGGRSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF",
"FUS_Mittal_Paper="=>"MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS",
"aSyn"=>"MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA",
"ht40"=>"MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL",
"A2"=>"GHMGRGGNFGFGDSRGGGGNFGPGPGSNFRGGSDGYGSGRGFGDGYNGYGGGPGGGNFGGSPGYGGGRGGYGGGGPGYGNQGGGYGGGYDNYGGGNYGSGNYNDFGNYNQQPSNYGPMKSGNFGGSRNMGGPYGGGNYGPGGSGGSGGYGGRSRY",
"A2D290V"=> "GHMGRGGNFGFGDSRGGGGNFGPGPGSNFRGGSDGYGSGRGFGDGYNGYGGGPGGGNFGGSPGYGGGRGGYGGGGPGYGNQGGGYGGGYDNYGGGNYGSGNYNVFGNYNQQPSNYGPMKSGNFGGSRNMGGPYGGGNYGPGGSGGSGGYGGRSRY",
"FUS12E"=>"MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA",
"p15PAF"=>"MVRTKADSVPGTYRKVVAARAPRKVLGSSTSATNSTSVSSRKAENKYAGGNPVCVRPTPKWQKGIGEFFRLSPKDSEKENQIPEEAGSSGLGKAKRKACPLQPDHTNDEKE",
"Sic1"=>"GSMTPSTPPRSRGTRYLAQPSGNTSSSALMQGQKTPQKPSQNLVPVTPSTTKSFKNAPLLAPPNSNMGMTSPFNGLTSPQRSPFPKSSVKR",
"Hst5"=>"DSHAKRHHGYKRKFHEKHHSHRGY",
"Hst52"=>"DSHAKRHHGYKRKFHEKHHSHRGYDSHAKRHHGYKRKFHEKHHSHRGY",
"ColNT"=>"MGSNGADNAHNNAFGGGKNPGIGNTSGAGSNGSASSNRGNSNGWSWSNKPHKNDGFHSDGSYHITFHGDNNSKPKPGGNSGNRGNNGDGASSHHHHHH",
"FhuA"=>"SESAWGPAATIAARQSATGTKTDTPIQKVPQSISVVTAEEMALHQPKSVKEALSYTPGVSVGTRGASNTYDHLIIRGFAAEGQSQNNYLNGLKLQGNFYNDAVIDPYMLERAEIMRGPVSVLYGKSSPGGLLNMVSKRPTTEPL",
"K25"=>"MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRL",
"OPN"=>"MHQDHVDSQSQEHLQQTQNDLASLQQTHYSSEENADVPEQPDFPDVPSKSQETVDDDDDDDNDSNDTDESDEVFTDFPTEAPVAPFNRGDNAGRGDSVAYGFRAKAHVVKASKIRKAARKLIEDDATTEDGDSQPAGLWWPKESREQNSRELPQHQSVENDSRPKFDSREVDGGDSKASAGVDSRESQGSVPAVDASNQTLESAEDAEDRHSIENNEVTR",
"PNt"=>"DWNNQSIVKTGERQHGIHIQGSDPGGVRTASGTTIKVSGRQAQGILLENPAAELQFRNGSVTSSGQLSDDGIRRFLGTVTVKAGKLVADHATLANVGDTWDDDGIALYVAGEQAQASIADSTLQGAGGVQIERGANVTVQRSAIVDGGLHIGALQSLQPEDLPPSRVVLRDTNVTAVPASGAPAAVSVLGASELTLDGGHITGGRAAGVAAMQGAVVHLQRATIRRGDALAGGAVPGGAVPGGAVPGGFGPGGFGPVLDGWYGVDVSGSSVELAQSIVEAPELGAAIRVGRGARVTVPGGSLSAPHGNVIETGGARRFAPQAAPLSITLQAGAH",
"FUS"=>"MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS",
"ACTR"=>"GTQNRPLLRNSLDDLVGPPSNLEGQSDERALLDQLHTLLSNTDATGLEEIDRALGIPELVNQGQALEPKQD",
"A1"=>"GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQSSNFGPMKGGNFGGRSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF",
"DmOrc1_expPh"=>"PLEIHLEQPEDNARPTRSSRKSLTAHRESKRSISARHDDTAGNKG##VEKRRRASMAASSSVEFIDVNSFICENKV#PIKIVGGRSVVRLSEKKNAPEINANYLPA#PL&EKNAKVE&PKSRASAARRNLNLSLDRGADTTADSDCLNY#IVQQ&PDPK&P#NDMKIKLRLSERRRSVRLA#MDVDPLSLEEAVQEPNAQGRKRLGVANGDIYH&P&KKSKEPLE#AAA&EQ&P#&RRK#ILK#ATSRLAEGTPRRSIHL#NIVEQRVFEDDEII#&PKRGRSKKTVQDNDEDY#PKK#VQKTPTRTRRSSTTTKTAT&PSKGIT&ATA&PM&PSQKMKKIRAGEL#PSMQQRTDLPAKDSSK",
"DmOrc1_thPh"=>"PLEIHLEQPEDNARPTRSSRKSLTAHRESKRSISARHDDTAGNKG##VEKRRRASMAASSSVEFIDVNSFICENKV#PIKIVGGRSVVRLSEKKNAPEINANYLPA#PL&EKNAKVE&PKSRASAARRNLNLSLDRGADTTADSDCLNY#IVQQ&PDPK&P#NDMKIKLRLSERRRSVRLA#MDVDPLSLEEAVQEPNAQGRKRLGVANGDIYH&P&KKSKEPLE#AAA&EQ&P#&RRK#ILK#ATSRLAEGTPRRSIHL#NIVEQRVFEDDEII#&PKRGRSKKTVQDNDEDY#PKK#VQKTPTRTRRSSTTTKTAT&PSKGIT&ATA&PM&PSQKMKKIRAGEL#PSMQQRTDLPAKDSSK",
"DmOrc1_noPh"=>"PLEIHLEQPEDNARPTRSSRKSLTAHRESKRSISARHDDTAGNKGSSVEKRRRASMAASSSVEFIDVNSFICENKVSPIKIVGGRSVVRLSEKKNAPEINANYLPASPLTEKNAKVETPKSRASAARRNLNLSLDRGADTTADSDCLNYSIVQQTPDPKTPSNDMKIKLRLSERRRSVRLASMDVDPLSLEEAVQEPNAQGRKRLGVANGDIYHTPTKKSKEPLESAAATEQTPSTRRKSILKSATSRLAEGTPRRSIHLSNIVEQRVFEDDEIISTPKRGRSKKTVQDNDEDYSPKKSVQKTPTRTRRSSTTTKTATTPSKGITTATATPMTPSQKMKKIRAGELSPSMQQRTDLPAKDSSK",
"htau40"=>"MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL",
"CorNID"=>"GPHMQVPRTHRLITLADHICQIITQDFARNQVPSQASTSTFQTSPSALSSTPVRTKTSSRYSPESQSQTVLHPRPGPRVSPENLVDKSRGSRPGKSPERSHIPSEPYEPISPPQGPAVHEKQDSMLLLSQRGVDPAEQRSDSRSPGSISYLPSFFTKLESTSPMVKSKKQEIFRKLNSSGGGDSDMAAAQPGTEIFNLPAVTTSGAVSSRSHSFADPASNLGLEDIIRKALMGSFDDKVEDHGVVMSHPVGIMPGSASTSVVTSSEARRDE",
"All-F"=>"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
"All-R"=>"RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR",
"All-P"=>"PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"
)

#uniprot + 
#https://phosphat.uni-hohenheim.de/phosphat.html
PhosSites = Dict{String, Vector{Int}}()
PhosSites["RS31"]  = [68,130,158,181,183,191,194,199,215,232,239,247,253,255,257]  #AT3G61860.1
PhosSites["RS31_Box"]  = [191,194,199,215]  #everything between 191:215
PhosSites["RS31_OUTSIDEBOX"]  = [68,158,181,183,232,239,247,253,255,257]   #everything outside 191:215

PhosSites["RS31a"] = [183,185,193,201,218,235,241,243]  #https://www.uniprot.org/uniprotkb/Q9ZPX8/feature-viewer
PhosSites["RS31a_Box"] = [193,201]   #everything between 193:201
PhosSites["RS31a_OUTSIDEBOX"] = [183,185,218,235,241,243]    #everything outside 193:201


PhosSites["RS40"]  = [89,151,154,178,193,195,203,211,216,221,222,230,232,241,262,267,269,278,282,288,298,306,308,309,317,320,326,335,340,342,345] #AT4G25500.1
PhosSites["RS41"]  = [171,176,177, 192,194, 202,205, 210, 219,220,229,231,239,254,272,274, 278,284,290,296,298,302, 309, 311, 313, 315,316,324,342,347,348,350,352] #https://www.psb.ugent.be/webtools/ptm-viewer/protein.php?id=AT5G52040.2


PhosDict=Dict('S'=>'#','T'=>'&','Y'=>'@')
function Phosphorylate(Sequence::String, Sites::Vector{Int})
    tmp = Vector{Char}(deepcopy(Sequence))
    for i in Sites
        tmp[i]=PhosDict[tmp[i]]
    end
    return join(tmp)
end

for name in ["RS31", "RS31a", "RS40", "RS41"]
    NameToSeq["$(name)_PHOSALL"] = Phosphorylate(NameToSeq[name], PhosSites[name])
end

NameToSeq["RS31_PHOS_BOX"] = Phosphorylate(NameToSeq["RS31"], PhosSites["RS31_Box"])
NameToSeq["RS31_PHOS_OUTSIDEBOX"] = Phosphorylate(NameToSeq["RS31"], PhosSites["RS31_OUTSIDEBOX"])
NameToSeq["RS31a_PHOS_BOX"] = Phosphorylate(NameToSeq["RS31a"], PhosSites["RS31a_Box"])
NameToSeq["RS31a_PHOS_OUTSIDEBOX"] = Phosphorylate(NameToSeq["RS31a"], PhosSites["RS31a_OUTSIDEBOX"])

### historical error at 157?
NameToSeq["RS41_PHOSALL_IDR2"] = NameToSeq["RS41_PHOSALL"][157:357]### according to Alpha-Fold 157-356 on mobiDB
NameToSeq["RS41_IDR2"] = NameToSeq["RS41"][157:357]### according to Alpha-Fold 157-356 on mobiDB
NameToSeq["RS41_PHOS6_IDR2"] = NameToSeq["RS41_PHOS6"][157:357]
NameToSeq["RS41_RRM"] = NameToSeq["RS41"][1:167] ### RS41_PHOS6_RRM is the same as RS41_RRM


NameToSeq["RS31a_PHOSALL_IDR2"] = NameToSeq["RS31a_PHOSALL"][167:end] ### according to Alpha-Fold 167-356 on mobiDB
NameToSeq["RS31_PHOSALL_IDR2"] = NameToSeq["RS31_PHOSALL"][167:end] ### according to Alpha-Fold 167-356 on mobiDB
NameToSeq["RS31_IDR2"] = NameToSeq["RS31"][167:length(NameToSeq["RS31a"])] ### according to Alpha-Fold 167-356 on mobiDB
NameToSeq["RS31a_IDR2"] = NameToSeq["RS31a"][167:length(NameToSeq["RS31a"])] ### according to Alpha-Fold 167-356 on mobiDB

NameToSeq["RS31_PHOS_BOX_IDR2"] = NameToSeq["RS31_PHOS_BOX"][167:end]
NameToSeq["RS31_PHOS_OUTSIDEBOX_IDR2"] = NameToSeq["RS31_PHOS_OUTSIDEBOX"][167:end]
NameToSeq["RS31a_PHOS_BOX_IDR2"] = NameToSeq["RS31a_PHOS_BOX"][167:end]
NameToSeq["RS31a_PHOS_OUTSIDEBOX_IDR2"] = NameToSeq["RS31a_PHOS_OUTSIDEBOX"][167:end]


NameToSeq["RS41_168_182"] = NameToSeq["RS41"][168:182] 


### From https://pubs.acs.org/doi/full/10.1021/bi800900d used for the Lindorf-Larsen R_Gs
### Domain Conformation of Tau Protein Studied by Solution Small-Angle X-ray Scattering 
NameToSeq["K32"] = "M"*NameToSeq["htau40"][198:394] #
NameToSeq["K16"] = "M"*NameToSeq["htau40"][198:372]
NameToSeq["K18"] = "M"*NameToSeq["htau40"][244:372]
NameToSeq["K27"] = "M"*NameToSeq["htau40"][198:276]*NameToSeq["htau40"][308:394]
NameToSeq["K17"] = "M"*NameToSeq["htau40"][198:276]*NameToSeq["htau40"][308:372]
NameToSeq["K19"] = "M"*NameToSeq["htau40"][244:276]*NameToSeq["htau40"][308:372]
NameToSeq["K44"] = "M"*NameToSeq["htau40"][1:44]*NameToSeq["htau40"][103:274]*NameToSeq["htau40"][308:372]
NameToSeq["K10"] = "M"*NameToSeq["htau40"][244:276]*NameToSeq["htau40"][308:length(NameToSeq["htau40"])] #
NameToSeq["K25"] = NameToSeq["htau40"][1:44]*NameToSeq["htau40"][103:243] #
NameToSeq["K23"] = NameToSeq["htau40"][1:44]*NameToSeq["htau40"][103:243]  *NameToSeq["htau40"][373:length(NameToSeq["htau40"])] 
NameToSeq["htau23"] = NameToSeq["htau40"][:44]*NameToSeq["htau40"][103:274]*NameToSeq["htau40"][308:length(NameToSeq["htau40"])]

#### Taken from file Modelling suggestions Stephan
NameToSeq["RS41_dCon1"] = "SRGNGYSPESPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon2"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSPDYGLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon3"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYVASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon5"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVNGSGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon4"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon12"] = "SRGNGYSPESPSPYRRERGSPDYGRGASPVAHKRERTSPDYGLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon13"] = "SRGNGYSPESPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon15"] = "SRGNGYSPESPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVNGSGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon125"] = "SRGNGYSPESPSPYRRERGSPDYGRGASPVAHKRERTSPDYGLSPDYKRDDRRRERVASPENGAVNGSGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon135"] = "SRGNGYSPESPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYVASPENGAVNGSGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dCon1235"] = "SRGNGYSPESPSPYRRERGSPDYGRGASPVAHKRERTSPDYGASPENGAVNGSGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dDom12"] = "SRGNGYSTSRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dTTPM"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRGPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dDom1"] = "SRGNGYSTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_dDom2"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_T0"] = "SRGNGYSPE"
NameToSeq["RS41_T1"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTS"
NameToSeq["RS41_T2"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRS"
NameToSeq["RS41_T3a"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSR"
NameToSeq["RS41_T4"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGSPDYGRGASPVAHKRERTSPDYGRGRRSPSPYKRARLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEV"

NameToSeq["RS41_B1"] = "SRGNGYSPERRRDRSPDRRRRSPSPYRRERGEPDYGRGAEPVAHKRERTEPDYGRGRRSPSPYKRARLEPDYKRDDRRRERVAEPENGAVRNREPRKGRGESRSPPPYEKRRESRSPPPYEKRRESRSPPPYEKRRERSRSRSKSSPENGQVESPGQIMEVEAGRGYDGADSPIRESPSRSPPAEE"
NameToSeq["RS41_B2"] = "ERGNGEEPERRRDRSPDRRRREPSPYRRERGSPDEGRGASPVAHKRERESPDYGRGRREPEPYKRARLSPDYKRDDRRRERVASPENGAVRNRSPRKGRGEEREPPPEEKRREEREPPPEEKRREEREPPPEEKRRERSREREKEEPENGQVEEPGQIMEVEAGRGYDGADEPIREEPEREPPAEE"
NameToSeq["RS41_B3"] = "ERGNGEEPERRRDRSPDRRRREPSPYRRERGEPDEGRGAEPVAHKREREEPDYGRGRREPEPYKRARLEPDYKRDDRRRERVAEPENGAVRNREPRKGRGEEREPPPEEKRREEREPPPEEKRREEREPPPEEKRRERSREREKEEPENGQVEEPGQIMEVEAGRGYDGADEPIREEPEREPPAEE"


### Taken from "Vorschlaege_Mutaten_RS31_RS31a
NameToSeq["RS31a_RS_11"] = "MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWAKDFQGERGKPRDGKAVSNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYASPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG"
NameToSeq["RS31a_RS_12"]="SPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG"
NameToSeq["RS31_RRM1"] = "MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWA"
NameToSeq["RS31_RRM2"] = "SNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYA"
NameToSeq["RS31_RRM_1"] = "KGERGRPRGDAKAP"
NameToSeq["RS31_RRM_2"] = "MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWASNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYA"
NameToSeq["RS31a_RRM1"] = "MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWA"
NameToSeq["RS31a_RRM2"] = "SNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYA"
NameToSeq["RS31a_RRM_IDR"] = "KDFQGERGKPRDGKAV"
NameToSeq["RS31a_RRM_2"] = "MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWASNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYA"
NameToSeq["RS31_RS"] = "LKDDDERDDRNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRARSPEYDRYKGPAAYERRRSPDYGRRSSDYGRQRSPGYDRYRSRSPVPRGRP"
NameToSeq["RS31a_RS"] = "LREAGEREDRYAGSRRRRSPSPVYRRRPSPDYTRRRSPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG"
NameToSeq["RS31_RS_1"] = "AAAAAAAAAAARNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRARSPEYDRYKGPAAYERRRSPDYGRRSSDYGRQRSPGYDRYRSRSPVPRGRP"
NameToSeq["RS31_RS_2"] = "RNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRARSPEYDRYKGPAAYERRRSPDYGRRSSDYGRQRSPGYDRYRSRSPVPRGRP"
NameToSeq["RS31_RS_3"] = "LKDDDERDDRNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRARSPEYDRYKGPAAYERRRSPDYGRRSSDY"
NameToSeq["RS31_RS_4"] = "LKDDDERDDRNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRAR"
NameToSeq["RS31_RS_5"] = "LKDDDERDD"
NameToSeq["RS31_RS_6"] = "MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWAKGERGRPRGDAKAPSNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYAAAAAAAAAAAARNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRARSPEYDRYKGPAAYERRRSPDYGRRSSDYGRQRSPGYDRYRSRSPVPRGRP"
NameToSeq["RS31_RS_7"] = "MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWAKGERGRPRGDAKAPSNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYARNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRARSPEYDRYKGPAAYERRRSPDYGRRSSDYGRQRSPGYDRYRSRSPVPRGRP"
NameToSeq["RS31_RS_8"] = "MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWAKGERGRPRGDAKAPSNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYALKDDDERDDRNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRARSPEYDRYKGPAAYERRRSPDYGRRSSDY"
NameToSeq["RS31_RS_9"] = "MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWAKGERGRPRGDAKAPSNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYALKDDDERDDRNGGRSPRRSLSPVYRRRPSPDYGRRPSPGQGRRPSPDYGRAR"
NameToSeq["RS31_RS_10"] = "MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWAKGERGRPRGDAKAPSNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYALKDDDERDD"
NameToSeq["RS31_RS_11"] = "MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWAKGERGRPRGDAKAPSNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYALKDDDERDDRNGGRSPRRSLSPVYRRRPSPDYGRRSPEYDRYKGPAAYERRRSPDYGRRSSDYGRQRSPGYDRYRSRSPVPRGRP"
NameToSeq["RS31_RS_12"] = "MRPVFVGNFEYETRQSDLERLFDKYGRVDRVDMKSGYAFVYFEDERDAEDAIRKLDNFPFGYEKRRLSVEWAKGERGRPRGDAKAPSNLKPTKTLFVINFDPIRTKEHDIEKHFEPYGKVTNVRIRRNFSFVQFETQEDATKALEATQRSKILDRVVSVEYASPEYDRYKGPAAYERRRSPDYGRRSSDYGRQRSPGYDRYRSRSPVPRGR"
NameToSeq["RS31_RS_13"] = "SPEYDRYKGPAAYERRRSPDYGRRSSDYGRQRSPGYDRYRSRSPVPRGR"

NameToSeq["RS31a_RS_1"] = "AAAAAAAAARYAGSRRRRSPSPVYRRRPSPDYTRRRSPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG"
NameToSeq["RS31a_RS_2"] = "RYAGSRRRRSPSPVYRRRPSPDYTRRRSPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG"
NameToSeq["RS31a_RS_3"] = "LREAGEREDRYAGSRRRRSPSPVYRRRPSPDYTRRRSPEYDRYKGPAPYERRKSPDYGRRSSDY"
NameToSeq["RS31a_RS_4"] = "LREAGEREDRYAGSRRRRSPSPVYRRRPSPDYTRRR"
NameToSeq["RS31a_RS_5"] = "LREAGERED"
NameToSeq["RS31a_RS_6"] = "MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWAKDFQGERGKPRDGKAVSNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYAAAAAAAAAARYAGSRRRRSPSPVYRRRPSPDYTRRRSPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG"
NameToSeq["RS31a_RS_7"] = "MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWAKDFQGERGKPRDGKAVSNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYARYAGSRRRRSPSPVYRRRPSPDYTRRRSPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG"
NameToSeq["RS31a_RS_8"] = "MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWAKDFQGERGKPRDGKAVSNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYALREAGEREDRYAGSRRRRSPSPVYRRRPSPDYTRRRSPEYDRYKGPAPYERRKSPDYGRRSSDY"
NameToSeq["RS31a_RS_9"] = "MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWAKDFQGERGKPRDGKAVSNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYALREAGEREDRYAGSRRRRSPSPVYRRRPSPDYTRRR"
NameToSeq["RS31a_RS_10"] = "MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWAKDFQGERGKPRDGKAVSNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYALREAGERED"
NameToSeq["RS31a_RS_11"] = "MRHVYVGNFDYDTRHSDLERLFSKFGRVKRVDMKSGYAFVYFEDERDAEDAIRRTDNTTFGYGRRKLSVEWAKDFQGERGKPRDGKAVSNQRPTKTLFVINFDPIRTRERDMERHFEPYGKVLNVRMRRNFAFVQFATQEDATKALDSTHNSKLLDKVVSVEYASPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG"
NameToSeq["RS31a_RS_12"] = "SPEYDRYKGPAPYERRKSPDYGRRSSDYGRARARSPGYDRSRSPIQRARG"

## MDP
NameToSeq["tia1"] = "GEDEMPKTLYVGNLSRDVTEALILQLFSQIGPCKNCKMIMDTAGNDPYCFVEFHEHRHAAAALAAMNGRKIMGKEVKVNWATTPSSQKLPQTGNHFHVFVGDLSPEITTEDIKAAFAPFGRISDARVVKDMATGKSKGYGFVSFFNKWDAENAIQQMGGQWLGGRQIRTNWATRKPPAPKSTYESNTKQLSYDEVVNQSSPSNCTVYCGGVTSGLTEQLMRQTFSPFGQIMEIRVFPDKGYSFVRFNSHESAAHAIVSVNGTTIEGHVVKCYWGK"
NameToSeq["gal3"] = "MADNFSLHDALSGSGNPNPQGWPGAWGNQPAGAGGYPGASYPGAYPGQAPPGAYPGQAPPGAYHGAPGAYPGAPAPGVYPGPPSGPGAYPSSGQPSAPGAYPATGPYGAPAGPLIVPYNLPLPGGVVPRMLITILGTVKPNANRIALDFQRGNDVAFHFNPRFNENNRRVIVCNTKLDNNWGREERQSVFPFESGKPFKIQVLVEPDHFKVAVNDAHLLQYNHRVKKLNEISKLGISGDIDLTSASYTMI"
NameToSeq["ubq2"] = "MASHHHHHHGAQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
NameToSeq["ubq3"] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
NameToSeq["ubq4"] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
end


