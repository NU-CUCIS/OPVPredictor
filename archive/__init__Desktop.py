# import flask, render_template to render the html pages, and request to be able to handle post requests between html pages
from flask import Flask, render_template, request
from sklearn.externals import joblib as pickle
# create a flask object
app = Flask(__name__)

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


# form is used for the id parameter in the html form
# this variable needs to be consistent accross html pages to pass the post parameters
htmlTitle = "Arindam Paul's"
template_text = "Testing"

#load things
model = loadData('bestModel')
dummyString = smiles = "Cc1csc(c2ccc(c3cc(C)c(c4cc5c(s4)c(c4ccc(C)cc4)c4ccsc4c5c4ccc(C)cc4)s3)c3nsnc23)c1"
smiles_form = 'form'




@app.route('/', methods = ['GET', 'POST'])
def homepage():

	return render_template("OPV.html", htmlTitle = htmlTitle, template_text = "", selectString = dummyString)

@app.route('/OPV.html', methods = ['GET', 'POST'])
def OPV():
	selectString=None
	print request.form.get('smiles', type=str)
	if selectString==None:
		selectString = dummyString
	else:
		print 'Not None'
	dummyFP = dummyAtomPairFromSmiles(selectString)
	output = str(model.predict(dummyFP)[0])



	return render_template("OPV.html", htmlTitle = htmlTitle, template_text = output, selectString = selectString)


if __name__ == "__main__":
	# to turn debug mode on use the following
	# app.run(port='6000')
	# app.run(debug=True,port=9043)
	from tornado.wsgi import WSGIContainer
	from tornado.httpserver import HTTPServer
	from tornado.ioloop import IOLoop
	http_server = HTTPServer(WSGIContainer(app))
	http_server.listen(9043)
	IOLoop.instance().start()
