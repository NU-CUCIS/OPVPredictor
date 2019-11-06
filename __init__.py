# import flask, render_template to render the html pages, and request to be able to handle post requests between html pages
from flask import Flask, render_template, request, Markup
import numpy,warnings
from sklearn.externals import joblib as pickle
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem,MACCSkeys,AtomPairs,Draw,SDMolSupplier,SDWriter
from rdkit.Chem.Draw import SimilarityMaps
from sklearn.ensemble import RandomForestRegressor
from matplotlib.colors import ColorConverter
from math import log10, floor
warnings.simplefilter("ignore")

# create a flask object
app = Flask(__name__)

def round_sig(x, sig=4):
    return round(x, sig-int(floor(log10(x)))-1)

def roundoff(x):
    return round(x,3)

def loadData(name,path="."):
    '''
    This loads a pickle file and returns the content which is a DICTIONARY object in our case.
    '''
    if ".pkl" in name:
            name = name.split(".pkl")[0]
    if "/" in name:
            name = name.split("/",1)[1]

    return pickle.load(path+"/"+name + '.pkl')

def dummyAtomPairFromSmiles(SMILE):
    fp = loadData('dummyFP')
    return fp


def getAtomPairFromSmiles(SMILE):
	mol = Chem.MolFromSmiles(SMILE)
	fp = SimilarityMaps.GetAPFingerprint(mol, fpType='bv')
	fpTranspose = numpy.transpose(fp)
	return fpTranspose

# form is used for the id parameter in the html form
# this variable needs to be consistent accross html pages to pass the post parameters
htmlTitle = "OPV Predictor"
template_text = "Testing"
#load things
model = loadData('maccsModel')

modelAtom = loadData('atomModel')

@app.route('/', methods = ['GET', 'POST'])
def homepage():
	return render_template("./index.html", htmlTitle = htmlTitle, template_text = "", selectString = "")

@app.route('/OPV.html', methods = ['GET', 'POST'])
def OPV():
    selectString = str(request.form.get('smiles', type=str))
    selectString2 = request.form.get('lumo', type=str)
    # selectString3 = request.form.get('filler', type=str)
    # selectString4 = request.form.get('current', type=str)
    # selectString5 = request.form.get('intensity', type=str)
    print 'In OPV'

    if selectString == "":
        selectString = "Invalid"
        print selectString1,"is empty"
    try:
        SMILE = selectString
        dummyFP = getAtomPairFromSmiles(SMILE)
        image = Draw.MolToImage(Chem.MolFromSmiles(SMILE), size=(250,250), fitImage=True, highlightAtoms=[1,2], highlightColor=ColorConverter().to_rgb('aqua'))
        image.save("static/molecule.png")
        print 'molecule saved'
        predictedAtom = model.predict(dummyFP)[0]
        # roundPredictedAtom = -round_sig(-predictedAtom)
        # print roundPredictedAtom
        B3LYP = roundoff(predictedAtom)
        PBE = roundoff(1.011*predictedAtom-0.007)
        BP86 = roundoff(0.906*predictedAtom-0.001)
        M06 = roundoff(1.066*predictedAtom-0.028)

        PBE_eV = roundoff((1.011*predictedAtom-0.007)*27.21)
        B3LYP_eV = roundoff(predictedAtom*27.21)
        BP86_eV = roundoff((0.906*predictedAtom-0.001)*27.21)
        M06_eV = roundoff((1.066*predictedAtom-0.028)*27.21)

        PBE_str = str(PBE)+" a.u. or "+str(PBE_eV)+" eV"
        B3LYP_str = str(B3LYP)+" a.u. or "+str(B3LYP_eV)+" eV"
        BP86_str = str(BP86)+" a.u. or "+str(BP86_eV)+" eV"
        M06_str = str(M06)+" a.u. or "+str(M06_eV)+" eV"

        # PBE_str = str(PBE)+" a.u. or "+str(PBE*27.21)+" eV"
        # B3LYP_str = str(B3LYP )+" a.u. or "+str(B3LYP *27.21)+" eV"
        # BP86_str = str(BP86)+" a.u. or "+str(BP86*27.21)+" eV"
        # M06_str = str(M06)+" a.u. or "+str(M06*27.21)+" eV"
        LUMO = abs(float(selectString2))
        HOMO = abs(predictedAtom)*27.21
        print "HOMO",HOMO
        # e = 1.60217662*pow(10,-19)
        # rec_e = 1/e
        # V = (energyDiff*pow(10,-20)- 0.3
        V_oc = HOMO - LUMO - 0.3
        print "Voltage", V_oc
        V_oc_string = str(round_sig(V_oc))+" V"
        print V_oc_string
        # FF = float(selectString3)
        # J = float(selectString4)
        # P_in = float(selectString5)
        # print FF, J, P_in
        # PCE_r = (V*FF*J)/P_in
        # PCE = 100*PCE_r
        # print PCE
        # PCE = round_sig(1/PCE_r)
        # pce = (-predicted)*100

	output=Markup("<table cellspacing='0' cellpading='0' style='width:90%'><tr><th>Functional</th><th>HOMO</th></tr><tr><td>PBE</td><td>"+PBE_str+"</td></tr><tr><td>B3LYP</td><td>"+B3LYP_str+"</td></tr><tr><td>BP86</td><td>"+BP86_str+"</td></tr><tr><td>M06</td><td>"+M06_str+"</td></tr></table>"+
    "<p>V_oc of donor-acceptor combination is "+V_oc_string+"</p>")
    #+"</td></tr><tr><td>PCE</td><td>"+str(PCE)
    #output = Markup(output1+output2)
    except:
        if selectString == "Invalid":
            selectString = ""
            output = Markup("<font color ='red'>You have not entered a SMILES formula</font>")

        else:
            output = Markup("<font color ='red'>Not valid SMILES Format</font>")

    return render_template("./OPV.html", htmlTitle = htmlTitle, template_text = output,  selectString = selectString)


if __name__ == "__main__":

	from tornado.wsgi import WSGIContainer
	from tornado.httpserver import HTTPServer
	from tornado.ioloop import IOLoop
	http_server = HTTPServer(WSGIContainer(app))
	http_server.listen(9043)
	IOLoop.instance().start()
