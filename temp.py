import matplotlib.pyplot as plt

r = 0.3
BR_Joe = [0,0,1-r,1,1]
Avon = [0,1-r,1-r,1-r,1]
BR_Avon = [0,1-r,1-r,1-r,1]
Joe = [0,0,1-r,1,1]
mark1_x = [-0.1,0]
mark1_y = [1-r,1-r]
mark2_x = [1-r,1-r]
mark2_y = [-0.1,0]
plt.plot(Joe,BR_Avon)
plt.plot(Avon,BR_Joe)
plt.legend([r'${\rm BR_{Avon}}(q)$',r'${\rm BR_{Joe}}(p)$'])
plt.scatter(1-r,1-r,marker='*',color='r')
plt.scatter(0,0,marker='o',color='k')
plt.scatter(1,1,marker='o',color='k')
plt.plot(mark1_x,mark1_y,'--',color='k')
plt.plot(mark2_x,mark2_y,'--',color='k')
plt.text(1-r+0.02,1-r+0.02,'msNE',color='r')
plt.text(0.02,0.02,'psNE',color='k')
plt.text(1-r-0.04,-0.17,r'$1-r$')
plt.text(-0.19,1-r-0.011,r'$1-r$')
plt.text(1-0.1,1-0.06,'psNE',color='k')
plt.xlabel(r'$p$')
plt.ylabel(r'$q$')
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])
plt.xticks([0,1],('Not Invest','Invest'))
plt.yticks([0,1],('Not Invest','Invest'))
plt.title(repr('$r='+str(r)+'$')[1:-1])
plt.show()