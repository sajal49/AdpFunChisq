from sklearn.gaussian_process import GaussianProcessRegressor
import numpy
import pandas

def gpr(x, y):

    # reshape
    x = numpy.asarray(x).reshape(-1, 1)
    y = numpy.asarray(y).reshape(-1, 1)
    
    # predict y
    gp = GaussianProcessRegressor().fit(x, y)
    y_predict = gp.predict(x)
    
    # create result
    x = numpy.array(x)
    ypr = numpy.array(y_predict)
    df = pandas.DataFrame(data={'x': x.flatten(), 'ypr': ypr.flatten()})
    
    return df
