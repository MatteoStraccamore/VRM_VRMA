import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata

class VRM:
    
    def __init__(
        self,
        Vreal,
        Iinput = 50.,
        phi0input = 0.,
        number_rings = 1,
        number_arch = 1,
        plot = False,
        save = False
    ):
        self.Vreal_los = Vreal[:,2]
        self.I = Iinput/180.*np.pi
        self.phi0 = phi0input/180.*np.pi
        self.passo = 1./number_rings
        self.passo2 = 2*np.pi/number_arch
        self.plot = plot
        self.save = save
        self.xp = Vreal[:,0]
        self.yp = Vreal[:,1]
        
    
    def _rescaled_R_theta(self,X,Y):
        r_projected=np.sqrt(X**2+Y**2)
        phi = np.nan_to_num(np.arccos(X/r_projected))
        phi = np.array([2*np.pi-phi[i] if Y[i] <= 0 else phi[i] for i in range(len(phi))])
        theta = np.arctan(np.tan(phi)/np.cos(self.I))
        theta = np.array([theta[i] + np.pi if ((phi[i]>np.pi/2) & (phi[i]<3.*np.pi/2))
                         else theta[i] + np.pi*2 if (phi[i]>3.*np.pi/2)
                         else 3.*np.pi/2 if (phi[i]==3.*np.pi/2)
                         else np.pi/2 if (phi[i]==np.pi/2)
                         else theta[i] for i in range(len(theta))])
        R = np.nan_to_num(r_projected*np.sqrt(np.cos(phi)**2+(np.sin(phi)**2)/np.cos(self.I)**2))
        R_rescaled = ((R-min(R))/(max(R)-min(R)))
        return R_rescaled, theta


    def _OFF(self):

        R_rescaled, theta = self._rescaled_R_theta(self.xp,self.yp)

        predictions = []
        a = np.sin(self.I)*np.cos(theta)
        b = np.sin(self.I)*np.sin(theta)
        ones = np.ones(shape=a.shape[0])
        X = np.vstack((ones,a,b)).T
        predictions.append(np.matmul(np.matmul(np.linalg.pinv(np.matmul(X.T,X)),X.T),self.Vreal_los))#,self.Vreal_los[array]))


        predictions = [[predictions[i][j] for i in range(len(predictions))] for j in range(3)]
        predictions = [item for sublist in predictions for item in sublist]
        predictions = np.array(predictions).reshape((1,len(predictions)))

        return predictions[0,0]
    
    
    #matrix; only for RTVR algorithm
    def _matrix(self):
        vsys = self._OFF()
        self.Vreal_los = self.Vreal_los - vsys
        print(vsys)
        R_rescaled, theta = self._rescaled_R_theta(self.xp,self.yp)
        theta_resc = theta

        predictions = []
        for i in np.arange(0, 1.0, self.passo):
            for ang in np.arange(0, 2.*np.pi, self.passo2):
                array = np.where((R_rescaled <=i+self.passo) & (R_rescaled >i) &
                                 (theta_resc <=ang+self.passo2) & (theta_resc >ang))

                a = np.sin(self.I)*np.cos(theta[array])
                b = np.sin(self.I)*np.sin(theta[array])
                X = np.vstack((a,b)).T
                predictions.append(np.matmul(np.matmul(np.linalg.pinv(np.matmul(X.T,X)),X.T),self.Vreal_los[array]))

        predictions = [[predictions[i][j] for i in range(len(predictions))] for j in range(2)]
        predictions = [item for sublist in predictions for item in sublist]
        predictions = np.array(predictions).reshape((1,len(predictions)))

        return predictions
    

    def _calculation_VLOS(self, predictions, file='no_name'):
    
        R_rescaled, theta = self._rescaled_R_theta(self.xp,self.yp)
        theta_resc = theta

        vT = predictions[0,0:int((predictions.shape[1])/2)]
        vR = predictions[0,int((predictions.shape[1])/2):int(predictions.shape[1])]
        j = 0
        Vlos = np.zeros((self.Vreal_los.shape[0]))
        for i in np.arange(0, 1.0, self.passo):
            for ang in np.arange(0, 2.*np.pi, self.passo2):
                array = np.where((R_rescaled <=i+self.passo) & (R_rescaled >=i) &
                                 (theta_resc <=ang+self.passo2) & (theta_resc >=ang))
                Vlos[array] = vT[j]*np.sin(self.I)*np.cos(theta[array])+ vR[j]*np.sin(self.I)*np.sin(theta[array])

                j = j+1
                
        if self.plot:
            plt.figure(figsize=(18,5))

            #REAL MAP
            plt.subplot(1, 2, 1)
            plt.title('Real data')

            # define grid.
            xi = np.linspace(np.amin(X),np.amax(X),100)
            yi = np.linspace(np.amin(Y),np.amax(Y),100)
            # grid the data.
            zi = griddata((self.xp, self.yp), self.Vreal_los, (xi[None,:], yi[:,None]), method='cubic')
            # contour the gridded data
            CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
            CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
            plt.colorbar() # draw colorbar

            # GENERATED MAP
            plt.subplot(1, 2, 2)
            plt.title('Generated data')

            # define grid.
            xi = np.linspace(np.amin(X),np.amax(X),100)
            yi = np.linspace(np.amin(Y),np.amax(Y),100)
            zi = griddata((self.xp, self.yp), Vlos, (xi[None,:], yi[:,None]), method='cubic') # grid the data
            # contour the gridded data
            CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
            CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
            plt.colorbar() # draw colorbar
            plt.show

        if self.save:
            plt.savefig("Vlos_%s.pdf" %(file))
            with open('Vlos_%s.txt' %(file), 'w') as f:
                f.write("xp'\typ'\ttheta\tVlos_real\tVlos_generated\tvR\tvT\tanello\tarco\n")
                j=0
                anello = 0
                for i in np.arange(0, 1.0, self.passo):
                    arco = 0
                    for ang in np.arange(0, 2*np.pi, self.passo2):
                        array = np.where((R_rescaled <=i+self.passo) & (R_rescaled >=i) &
                                         (theta_resc <=ang+self.passo2) & (theta_resc >=ang))
                        for k in range(0,array[0].shape[0]):
                            f.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\t%i\n" %(X[array[0][k]],Y[array[0][k]],
                                                                             theta[array[0][k]],
                                                                             self.Vreal_los[array[0][k]],
                                                                             Vlos[array[0][k]],
                                                                             vR[j],vT[j],anello,arco))
                        j=j+1
                        arco = arco + 1
                    anello = anello+1
                    
        return Vlos
