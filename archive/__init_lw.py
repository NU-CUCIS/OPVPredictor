# import flask, render_template to render the html pages, and request to be able to handle post requests between html pages
from flask import Flask, render_template, request, Markup
import numpy,warnings
from sklearn.externals import joblib as pickle
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem,MACCSkeys,AtomPairs
from rdkit.Chem.Draw import SimilarityMaps
from sklearn.ensemble import RandomForestRegressor
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


@app.route('/', methods = ['GET', 'POST'])
def homepage():
	return render_template("./index.html", htmlTitle = htmlTitle, template_text = "", selectString = "")

@app.route('/OPV.html', methods = ['GET', 'POST'])
def OPV():
    selectString = str(request.form.get('smiles', type=str))
    #selectString2 = request.form.get('lumo', type=str)
    #selectString3 = request.form.get('filler', type=str)
    #selectString4 = request.form.get('current', type=str)
    #selectString5 = request.form.get('intensity', type=str)

    if selectString == "":
        selectString = "Invalid"
        # print selectString1,"is empty"
    try:
        dummyFP = getAtomPairFromSmiles(selectString)
        predicted = model.predict(dummyFP)[0]
        #print predicted
        roundPredicted = -round_sig(-predicted)
        print predicted, roundPredicted
        #lumo = float(selectString2)
        #energyDiff = predicted - lumo
        #rec_V = energyDiff - 0.3
        #e = 1.60217662*pow(10,-19)
        #rec_e = 1/e
        #V = rec_e*energyDiff*pow(10,-20)
        #roundV = round_sig(V)
        #print V,roundV
        #FF = float(selectString3)
        #J = float(selectString4)
        #P_in = float(selectString5)

        #PCE_r = (V*FF*J)/P_in
        #PCE = round_sig(1/PCE_r)
        #print PCE
        # pce = (-predicted)*100
        #output = Markup("The predicted HOMO value is <font color ='red'>"+str(roundPredicted)+"eV</font>, the voltage is <font color ='red'>"+str(roundV)+" V</font> and the predicted power conversion efficiency is <font color ='red'>"+str(PCE)+"%</font>")
	output = Markup("The predicted HOMO value is <font color ='red'>"+str(roundPredicted)+"eV</font>")
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
