"""Reproduce the corrected vector diagrams for the lecture section.

Requires NumPy and Matplotlib. Run from any directory; PDF outputs are written
to the media paths referenced by the LaTeX file. All points are schematic.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse, Arc

ROOT = Path(__file__).resolve().parent
plt.rcParams.update({'font.family': 'DejaVu Serif', 'font.size': 11,
                     'pdf.fonttype': 42, 'mathtext.fontset': 'dejavuserif'})
blue, red, green = '#006FA5', "#5d1007", "#275850"

def flow():
    fig, ax = plt.subplots(figsize=(12, 4.4))
    ax.set(xlim=(0, 12), ylim=(0, 4.4))
    ax.axis('off')
    boxes = [(0.1, 1.2, 3.0, 2.6), (4, 1.2, 4, 2.6), (8.9, 1.2, 3, 2.6)]
    for x, y, w, h in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h, boxstyle='round,pad=0.06,rounding_size=0.13',
                                   edgecolor=blue, facecolor='#f6f8fb', linewidth=1.5))
    ax.text(1.6,3.4,'Indirekte Messgrößen',ha='center',weight='bold',fontsize=13)
    ax.text(1.6,2.77,r'$Y_1,\ldots,Y_M$',ha='center',fontsize=15)
    ax.text(1.6,2.18,'Beschreibung durch',ha='center')
    ax.text(1.6,1.73,r'unbekannte Parameter $\mathbf{p}_0$',ha='center')
    ax.text(6,3.4,'Vorwärtsmodell',ha='center',weight='bold',fontsize=13)
    ax.text(6,2.85,r'$\mathbf{X}_j=\mathbf{f}_j(\mathbf{p}_0)+\mathbf{E}_j$',ha='center',fontsize=15)
    ax.text(6,2.27,r'Physik: $\mathbf{f}_j(\mathbf{p}_0)$',ha='center')
    ax.text(6,1.72,r'Statistik: Verteilung der $\mathbf{E}_j$',ha='center')
    ax.text(10.4,3.4,'Direkte Messgrößen',ha='center',weight='bold',fontsize=13)
    ax.text(10.4,2.77,r'$X_{1,j},\ldots,X_{N,j}$',ha='center',fontsize=15)
    ax.text(10.4,2.18,'Realisierte Messwerte',ha='center')
    ax.text(10.4,1.73,r'$\mathbf{x}_1,\ldots,\mathbf{x}_J$',ha='center',fontsize=14)
    for a,b in [(3.2,3.88),(8.12,8.78)]:
        ax.add_patch(FancyArrowPatch((a,2.55),(b,2.55),arrowstyle='-|>',mutation_scale=20,color=blue,lw=2))
    ax.plot([10.4,10.4],[1.12,.64],color=green,lw=1.7)
    ax.plot([1.6,1.6],[.64,1.05],color=green,lw=1.7)
    ax.add_patch(FancyArrowPatch((10.4,.64),(1.6,.64),arrowstyle='-|>',mutation_scale=18,color=green,lw=1.7))
    ax.text(6,.17,r'Inverses Problem: Parameterschätzung $\widehat{\mathbf{p}}$ aus den Messwerten',
            ha='center',color=green,fontsize=12)
    fig.savefig(ROOT/'04_vorlesung/media/Modell1.pdf',bbox_inches='tight',pad_inches=.04)
    plt.close(fig)

def projections():
    fig, axes = plt.subplots(1,2,figsize=(11.8,6.0))
    slope, intercept = .55, .50
    alpha=np.arctan(slope)
    n=np.array([-np.sin(alpha),np.cos(alpha)])
    d=intercept*np.cos(alpha)
    x=np.array([2.,3.6])
    eps=n@x-d
    zp=x-eps*n
    for i,(ax,S) in enumerate(zip(axes,[np.eye(2),np.diag([1.,.25])])):
        r=eps*(S@n)/(n@S@n)
        z=x-r
        q=eps**2/(n@S@n)
        grid=np.array([-.55,5.35])
        ax.plot(grid,slope*grid+intercept,color=blue,lw=2.2,zorder=1)
        ax.add_patch(Ellipse(x,2*np.sqrt(q*S[0,0]),2*np.sqrt(q*S[1,1]),
                            facecolor=green,alpha=.07,edgecolor='none'))
        ax.add_patch(Ellipse(x,2*np.sqrt(q*S[0,0]),2*np.sqrt(q*S[1,1]),
                            facecolor='none',edgecolor=green,lw=1.5))
        ax.scatter(*x,s=75,facecolor='#efb93e',edgecolor=blue,zorder=5)
        ax.text(*(x+[-.12,.19]),r'$\mathbf{x}_j$',ha='right',fontsize=14)
        if i:
            ax.plot([x[0],zp[0]],[x[1],zp[1]],'--',color='#777777',lw=1.5)
            ax.scatter(*zp,s=38,facecolor='white',edgecolor='#777777',zorder=5)
            ax.text(*(zp+[-.17,-.32]),r'$\mathbf{z}_{\perp,j}$',fontsize=12,ha='right')
            ax.text(2.28,2.57,r'$\varepsilon_j$',fontsize=13,color='#666666')
        ax.annotate('',xy=x,xytext=z,arrowprops={'arrowstyle':'->','lw':2.1,'color':red})
        ax.scatter(*z,s=42,facecolor=red,edgecolor='white',zorder=6)
        if i:
            ax.text(*(z+[.08,-.20]),r'$\mathbf{z}_j^*$',fontsize=13)
            ax.text(3.21,3.42,r'$\boldsymbol{\varepsilon}_j^*$',fontsize=14,color=red)
        else:
            ax.text(*(z+[.14,-.10]),r'$\mathbf{z}_{\perp,j}$',fontsize=13)
            ax.text(2.65,3.06,r'$\boldsymbol{\varepsilon}_{\perp,j}=\varepsilon_j\mathbf{n}$',fontsize=12,color=red)
            tang=np.array([np.cos(alpha),np.sin(alpha)])
            pts=np.array([zp+.18*n,zp+.18*n+.18*tang,zp+.18*tang])
            ax.plot(pts[:,0],pts[:,1],color=red,lw=1)
        ax.set(xlim=(-1.0,5.4),ylim=(-.25,5.65),aspect='equal')
        ax.spines[['top','right']].set_visible(False)
        ax.spines[['left','bottom']].set_position('zero')
        ax.spines[['left','bottom']].set_color('#777777')
        ax.set_xticks([]);ax.set_yticks([])
        ax.text(5.20,-.17,r'$x_1$',fontsize=13)
        ax.text(-.34,5.42,r'$x_2$',fontsize=13)
        ax.text(4.16,2.18,r'$\mathbf{n}^{\mathsf{T}}\mathbf{z}=d$',color=blue,fontsize=11,rotation=np.degrees(alpha))
        ax.set_title(['Isotrope Koordinatenkovarianz','Unterschiedliche Koordinatenvarianzen'][i],fontsize=12,pad=16)
        ax.text(.5,-.085,[r'$\Sigma=\sigma^2 I$   (hier $\sigma=1$)',
                          r'$\Sigma=\operatorname{diag}(1,\,0{,}25)$'][i],
                transform=ax.transAxes,ha='center',fontsize=12)
    fig.subplots_adjust(left=.04,right=.985,bottom=.14,top=.91,wspace=.16)
    fig.savefig(ROOT/'04_vorlesung/media/ResiduenVarianzen.pdf',bbox_inches='tight',pad_inches=.08)
    plt.close(fig)

if __name__=='__main__':
#    (ROOT/'04_vorlesung/media').mkdir(parents=True, exist_ok=True)
    flow()
    projections()
