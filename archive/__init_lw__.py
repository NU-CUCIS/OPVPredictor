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
model = loadData('bestModel')

modelAtom = loadData('atomModel')

@app.route('/', methods = ['GET', 'POST'])
def homepage():
	return render_template("./index.html", htmlTitle = htmlTitle, template_text = "", selectString = "")

@app.route('/OPV.html', methods = ['GET', 'POST'])
def OPV():
    selectString = str(request.form.get('smiles', type=str))
    print 'In OPV'

    if selectString == "":
        selectString = "Invalid"
        print selectString1,"is empty"
    try:
        SMILE = selectString
        dummyFP = getAtomPairFromSmiles(SMILE)
        image = Draw.MolToImage(Chem.MolFromSmiles(SMILE),size=(1000, 1000), highlightAtoms=[1,2], highlightColor=ColorConverter().to_rgb('aqua'))
        image.save("static/molecule.png")
        print 'molecule saved'
        predictedAtom = model.predict(dummyFP)[0]
        roundPredictedAtom = -round_sig(-predictedAtom)
        B3LYP = roundPredictedAtom
	PBE = 1.011*roundPredictedAtom-0.007
        BP86 = 0.906*roundPredictedAtom-0.001
        M06 = 1.066*roundPredictedAtom-0.028
	output=Markup("<table style='width:100%'><tr><th>Functional</th><th>HOMO (in eV)</th></tr><tr><td>PBE</td><td>"+str(PBE)+"</td></tr><tr><td>B3LYP</td><td>"+str(B3LYP)+"</td></tr><tr><td>BP86</td><td>"+str(BP86)+"</td></tr><tr><td>M06</td><td>"+str(M06)+"</td></tr></table>")

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
