/* ===================== MODULO 1 - SISTEMAS DE ECUACIONES LINEALES =====================
   Metodos implementados desde cero: Jacobi, Gauss-Seidel, SOR, LU, Gradiente Conjugado */
const PLANTAS=[{name:"Planta Achacachi",code:"x\u2081 \u00b7 zona norte"},{name:"Planta Senkata",code:"x\u2082 \u00b7 El Alto"},{name:"Planta Patacamaya",code:"x\u2083 \u00b7 zona sur"}];
const HOSP=["H. Municipal Los Andes","Hospital Holand\u00e9s","Hospital La Paz"];
const DEFAULT_A=[[0.8,0.1,0.1],[0.1,0.7,0.2],[0.1,0.2,0.6]];
const DEFAULT_B=[120,150,180];
const M1_METHODS={
  gauss:{label:'Gauss-Seidel',kind:'iter',desc:'Iterativo. Reutiliza los valores ya actualizados en el mismo paso. Converge si la matriz es diagonalmente dominante.'},
  jacobi:{label:'Jacobi',kind:'iter',desc:'Iterativo. Calcula todas las inc\u00f3gnitas con los valores del paso anterior; suele ser m\u00e1s lento que Gauss-Seidel. Converge si la matriz es diagonalmente dominante.'},
  sor:{label:'SOR',kind:'iter',desc:'Iterativo. Es Gauss-Seidel con un factor \u03c9 que acelera o frena la convergencia. Requiere 0 < \u03c9 < 2.'},
  lu:{label:'Descomposici\u00f3n LU',kind:'direct',desc:'Directo. Factoriza A = L\u00b7U y resuelve en un solo paso. Funciona con cualquier matriz no singular (determinante \u2260 0).'},
  cg:{label:'Gradiente Conjugado',kind:'iter',desc:'Iterativo y muy eficiente, pero solo funciona si la matriz es sim\u00e9trica y definida positiva.'}
};
let convChart=null;

function buildInputs(){
  const mg=document.getElementById('matrixGrid'),bg=document.getElementById('bvecGrid'); mg.innerHTML='';bg.innerHTML='';
  for(let i=0;i<3;i++){
    for(let j=0;j<3;j++){const e=document.createElement('input');e.type='number';e.step='0.01';e.className='num-in';e.id='a_'+i+'_'+j;e.value=DEFAULT_A[i][j];mg.appendChild(e);}
    const b=document.createElement('input');b.type='number';b.step='1';b.className='num-in';b.id='b_'+i;b.value=DEFAULT_B[i];b.title=HOSP[i];bg.appendChild(b);
  }
}
function readSystem(){const A=[],b=[];for(let i=0;i<3;i++){A.push([]);for(let j=0;j<3;j++)A[i].push(parseFloat(document.getElementById('a_'+i+'_'+j).value)||0);b.push(parseFloat(document.getElementById('b_'+i).value)||0);}return{A,b};}

// ---- utilidades de matriz ----
function isDominant(A){let dom=true,strict=false;for(let i=0;i<3;i++){const d=Math.abs(A[i][i]);let off=0;for(let j=0;j<3;j++)if(j!==i)off+=Math.abs(A[i][j]);if(d<off)dom=false;if(d>off)strict=true;}return dom&&strict;}
function isSymmetric(A){for(let i=0;i<3;i++)for(let j=i+1;j<3;j++)if(Math.abs(A[i][j]-A[j][i])>1e-9)return false;return true;}
function det3(A){return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);}
function isPosDef(A){if(!isSymmetric(A))return false;const m1=A[0][0],m2=A[0][0]*A[1][1]-A[0][1]*A[1][0],m3=det3(A);return m1>0&&m2>0&&m3>0;}
function dot(a,b){let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s;}
function matVec(A,v){return A.map(r=>r.reduce((s,a,j)=>s+a*v[j],0));}
function relErr(xn,xo){let e=0;for(let i=0;i<xn.length;i++){const d=Math.abs(xn[i])>1e-12?Math.abs(xn[i]):1;e=Math.max(e,Math.abs(xn[i]-xo[i])/d*100);}return e;}

// ---- METODOS ----
function mJacobi(A,b,tol,maxIter){
  let x=[0,0,0];const rows=[];let k=0,err=Infinity,div=false;
  while(k<maxIter){const xn=[0,0,0];
    for(let i=0;i<3;i++){let s=b[i];for(let j=0;j<3;j++)if(j!==i)s-=A[i][j]*x[j];xn[i]=s/A[i][i];}
    err=relErr(xn,x);x=xn;k++;rows.push({k,x:[...x],err});
    if(!x.every(Number.isFinite)||err>1e10){div=true;break;} if(err<tol)break;}
  return {type:'iter',x,rows,iterations:k,finalErr:err,converged:err<tol&&!div,diverged:div};
}
function mGaussSeidel(A,b,tol,maxIter){
  let x=[0,0,0];const rows=[];let k=0,err=Infinity,div=false;
  while(k<maxIter){const xo=[...x];
    for(let i=0;i<3;i++){let s=b[i];for(let j=0;j<3;j++)if(j!==i)s-=A[i][j]*x[j];x[i]=s/A[i][i];}
    err=relErr(x,xo);k++;rows.push({k,x:[...x],err});
    if(!x.every(Number.isFinite)||err>1e10){div=true;break;} if(err<tol)break;}
  return {type:'iter',x,rows,iterations:k,finalErr:err,converged:err<tol&&!div,diverged:div};
}
function mSOR(A,b,w,tol,maxIter){
  let x=[0,0,0];const rows=[];let k=0,err=Infinity,div=false;
  while(k<maxIter){const xo=[...x];
    for(let i=0;i<3;i++){let s=b[i];for(let j=0;j<3;j++)if(j!==i)s-=A[i][j]*x[j];const gs=s/A[i][i];x[i]=(1-w)*x[i]+w*gs;}
    err=relErr(x,xo);k++;rows.push({k,x:[...x],err});
    if(!x.every(Number.isFinite)||err>1e10){div=true;break;} if(err<tol)break;}
  return {type:'iter',x,rows,iterations:k,finalErr:err,converged:err<tol&&!div,diverged:div};
}
function mLU(A,b){
  const n=3,U=A.map(r=>r.slice()),L=[[1,0,0],[0,1,0],[0,0,1]],pb=b.slice();
  for(let k=0;k<n;k++){
    let p=k,mx=Math.abs(U[k][k]);
    for(let i=k+1;i<n;i++)if(Math.abs(U[i][k])>mx){mx=Math.abs(U[i][k]);p=i;}
    if(mx<1e-12)return {type:'lu',singular:true};
    if(p!==k){const t=U[k];U[k]=U[p];U[p]=t;const tb=pb[k];pb[k]=pb[p];pb[p]=tb;for(let j=0;j<k;j++){const tl=L[k][j];L[k][j]=L[p][j];L[p][j]=tl;}}
    for(let i=k+1;i<n;i++){const f=U[i][k]/U[k][k];L[i][k]=f;for(let j=k;j<n;j++)U[i][j]-=f*U[k][j];}
  }
  const y=[0,0,0];for(let i=0;i<n;i++){let s=pb[i];for(let j=0;j<i;j++)s-=L[i][j]*y[j];y[i]=s;}
  const x=[0,0,0];for(let i=n-1;i>=0;i--){let s=y[i];for(let j=i+1;j<n;j++)s-=U[i][j]*x[j];x[i]=s/U[i][i];}
  return {type:'lu',singular:false,x,L,U};
}
function mCG(A,b,tol,maxIter){
  let x=[0,0,0],r=b.slice(),p=r.slice(),rs=dot(r,r);const rows=[];let k=0;const nb=Math.sqrt(dot(b,b))||1;
  while(k<maxIter){
    const Ap=matVec(A,p),alpha=rs/dot(p,Ap);
    for(let i=0;i<3;i++){x[i]+=alpha*p[i];r[i]-=alpha*Ap[i];}
    const rsn=dot(r,r);k++;const res=Math.sqrt(rsn)/nb*100;rows.push({k,x:[...x],err:res});
    if(Math.sqrt(rsn)<1e-12||res<tol){rs=rsn;break;}
    const beta=rsn/rs;for(let i=0;i<3;i++)p[i]=r[i]+beta*p[i];rs=rsn;
  }
  return {type:'iter',cg:true,x,rows,iterations:k,finalErr:rows.length?rows[rows.length-1].err:0,converged:true,diverged:false};
}

// ---- badges de requisito ----
function bOk(t){return '<span class="dom-badge dom-ok"><span class="dot" style="background:var(--teal)"></span> '+t+'</span>';}
function bBad(t){return '<span class="dom-badge dom-bad"><span class="dot" style="background:var(--crisis)"></span> '+t+'</span>';}
function renderRequirement(method,A){
  const el=document.getElementById('domBadge');
  if(method==='cg'){const sym=isSymmetric(A);el.innerHTML=(sym&&isPosDef(A))?bOk('Sim\u00e9trica y definida positiva \u00b7 apta para Gradiente Conjugado'):bBad(sym?'Sim\u00e9trica pero NO definida positiva':'NO sim\u00e9trica \u00b7 Gradiente Conjugado no aplicable');}
  else if(method==='lu'){const d=det3(A);el.innerHTML=Math.abs(d)>1e-9?bOk('Determinante = '+d.toFixed(3)+' \u2260 0 \u00b7 LU aplicable'):bBad('Determinante \u2248 0 \u00b7 matriz singular, LU no aplicable');}
  else el.innerHTML=isDominant(A)?bOk('Diagonalmente dominante \u00b7 convergencia garantizada'):bBad('NO dominante \u00b7 puede no converger');
}

// ---- render ----
function hideChart(){document.getElementById('m1_chartPanel').style.display='none';}
function showChart(res,method){
  const panel=document.getElementById('m1_chartPanel');panel.style.display='';
  document.getElementById('m1_chartLabel').textContent=res.cg?'Convergencia del residuo':'Convergencia del error';
  const rows=res.rows,ctx=document.getElementById('convChart');
  if(convChart)convChart.destroy();
  convChart=new Chart(ctx,{type:'line',data:{labels:rows.map(r=>r.k),datasets:[{label:res.cg?'Residuo (%)':'Error relativo (%)',data:rows.map(r=>Math.max(r.err,1e-12)),borderColor:'#E4A64D',backgroundColor:'rgba(228,166,77,.12)',borderWidth:2,pointRadius:3,pointBackgroundColor:'#E4A64D',tension:.25,fill:true}]},
    options:{responsive:true,plugins:{legend:{labels:{color:'#B6AC94',font:{family:'IBM Plex Mono',size:11}}}},
      scales:{y:{type:'logarithmic',grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}}},x:{grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}}}}}});
}
function renderSolution(res){
  document.getElementById('solBox').innerHTML=PLANTAS.map((p,i)=>'<div class="x-row"><div class="name">'+p.name+'<small>'+p.code+'</small></div><div class="qty">'+res.x[i].toFixed(2)+' <span>cil.</span></div></div>').join('');
}
function renderIterTable(res){
  document.getElementById('m1_detailLabel').textContent='Tabla de iteraciones';
  const eh=res.cg?'residuo %':'error %';
  let h='<div class="tbl-scroll"><table class="tbl-iter"><thead><tr><th>k</th><th>x\u2081</th><th>x\u2082</th><th>x\u2083</th><th>'+eh+'</th></tr></thead><tbody>';
  h+=res.rows.map(r=>'<tr><td>'+r.k+'</td><td>'+r.x[0].toFixed(4)+'</td><td>'+r.x[1].toFixed(4)+'</td><td>'+r.x[2].toFixed(4)+'</td><td>'+r.err.toExponential(2)+'</td></tr>').join('');
  h+='</tbody></table></div>';
  document.getElementById('m1_detail').innerHTML=h;
}
function matHTML(M,title){let h='<div class="lu-mat"><div class="lu-title">'+title+'</div><table class="lu-table">';for(let i=0;i<3;i++){h+='<tr>';for(let j=0;j<3;j++)h+='<td>'+M[i][j].toFixed(3)+'</td>';h+='</tr>';}return h+'</table></div>';}
function renderLU(res){
  document.getElementById('m1_detailLabel').textContent='Descomposici\u00f3n A = L \u00b7 U';
  document.getElementById('m1_detail').innerHTML='<div class="lu-wrap">'+matHTML(res.L,'L (inferior)')+matHTML(res.U,'U (superior)')+'</div>';
}
function renderInterp(res,method,A){
  const box=document.getElementById('interpBox'),txt=document.getElementById('interpText');box.style.display='block';
  const name=M1_METHODS[method].label;let h='';
  if(res.type==='lu'){h='<b>'+name+'</b> es un m\u00e9todo <b>directo</b>: factoriza A = L\u00b7U y resuelve por sustituci\u00f3n, sin iteraciones. Da la soluci\u00f3n exacta (salvo redondeo) siempre que la matriz no sea singular. Soluci\u00f3n: x\u2081='+res.x[0].toFixed(2)+', x\u2082='+res.x[1].toFixed(2)+', x\u2083='+res.x[2].toFixed(2)+'.';}
  else if(res.converged){h='<b>'+name+'</b> convergi\u00f3 en <b>'+res.iterations+'</b> iteraciones (error final '+res.finalErr.toExponential(2)+(res.cg?' % de residuo). ':' %). ');h+=res.cg?'El Gradiente Conjugado aprovecha que la matriz es sim\u00e9trica y definida positiva; para n inc\u00f3gnitas converge en a lo sumo n pasos.':(isDominant(A)?'Como la matriz es diagonalmente dominante, la convergencia estaba garantizada.':'La matriz no es estrictamente dominante, pero aun as\u00ed convergi\u00f3.');}
  else if(res.diverged){h='<b>'+name+'</b> <b style="color:var(--crisis)">divergi\u00f3</b>: al no cumplir la matriz la condici\u00f3n de dominancia, las correcciones crecen sin control. Prueba con LU, que no necesita esa condici\u00f3n.';}
  else{h='<b>'+name+'</b> no alcanz\u00f3 la tolerancia en '+res.iterations+' iteraciones (error '+res.finalErr.toExponential(2)+' %). Converge muy lento: sube el m\u00e1ximo de iteraciones o prueba otro m\u00e9todo.';}
  txt.innerHTML=h;
}
function renderPlain(res,method){
  const el=document.getElementById('plainConclusion');if(!el)return;const name=M1_METHODS[method].label;
  if(res.type==='lu'){el.innerHTML='Usando <b>'+name+'</b> (m\u00e9todo directo), el sistema se resolvi\u00f3 en un solo paso: las plantas deber\u00edan despachar ~<b>'+Math.round(res.x[0])+'</b>, <b>'+Math.round(res.x[1])+'</b> y <b>'+Math.round(res.x[2])+'</b> cilindros (Achacachi, Senkata, Patacamaya).';return;}
  if(res.converged){el.innerHTML='Con <b>'+name+'</b>, el sistema convergi\u00f3 en <b>'+res.iterations+'</b> '+(res.cg?'iteraciones':'ajustes')+': las plantas deber\u00edan despachar ~<b>'+Math.round(res.x[0])+'</b>, <b>'+Math.round(res.x[1])+'</b> y <b>'+Math.round(res.x[2])+'</b> cilindros.';}
  else if(res.diverged){el.innerHTML='Con <b>'+name+'</b> el sistema <b>no converge</b>: los valores se disparan en vez de estabilizarse. La matriz no cumple las condiciones de este m\u00e9todo; prueba con <b>LU</b>.';}
  else{el.innerHTML='Con <b>'+name+'</b> no se alcanz\u00f3 la tolerancia en el m\u00e1ximo de iteraciones. Aumenta el m\u00e1ximo o prueba otro m\u00e9todo.';}
}
function renderBlocked(msg,method){
  document.getElementById('solBox').innerHTML='<div class="blocked-msg">\u26a0 No se puede resolver con '+M1_METHODS[method].label+'<p>'+msg+'</p></div>';
  document.getElementById('m1_detailLabel').textContent='Detalle';
  document.getElementById('m1_detail').innerHTML='<p class="dim mb-0" style="padding:.8rem">Sin resultados: el m\u00e9todo elegido no es aplicable a esta matriz.</p>';
  hideChart();
  const box=document.getElementById('interpBox');box.style.display='block';
  document.getElementById('interpText').innerHTML=msg+' Cada m\u00e9todo num\u00e9rico exige que la matriz cumpla ciertas condiciones; cuando no se cumplen, hay que elegir otro.';
  document.getElementById('plainConclusion').innerHTML='Con este m\u00e9todo no se puede resolver el sistema tal como est\u00e1, porque la matriz no cumple el requisito que el m\u00e9todo necesita. Elige otro m\u00e9todo en el men\u00fa.';
}

function solve(){
  const {A,b}=readSystem();
  const method=document.getElementById('m1_method').value;
  renderRequirement(method,A);
  const tol=parseFloat(document.getElementById('tol').value)||0.01;
  const mx=parseInt(document.getElementById('maxit').value)||100;
  if(method==='cg'){
    if(!isSymmetric(A)){renderBlocked('El Gradiente Conjugado solo funciona con matrices <b>sim\u00e9tricas</b> (a\u1d62\u2c7c = a\u2c7c\u1d62). Esta matriz no lo es, as\u00ed que no se puede aplicar. Prueba con <b>LU</b> o <b>Gauss-Seidel</b>.',method);return;}
    if(!isPosDef(A)){renderBlocked('La matriz es sim\u00e9trica pero <b>no definida positiva</b> (alg\u00fan menor principal es \u2264 0), condici\u00f3n que el Gradiente Conjugado necesita. Prueba con <b>LU</b>.',method);return;}
  }
  let res;
  if(method==='jacobi')res=mJacobi(A,b,tol,mx);
  else if(method==='gauss')res=mGaussSeidel(A,b,tol,mx);
  else if(method==='sor'){const w=parseFloat(document.getElementById('m1_omega').value);if(!(w>0&&w<2)){renderBlocked('SOR requiere un factor de relajaci\u00f3n <b>\u03c9 entre 0 y 2</b>. Ajusta ese valor.',method);return;}res=mSOR(A,b,w,tol,mx);}
  else if(method==='lu'){res=mLU(A,b);if(res.singular){renderBlocked('La matriz es <b>singular</b> (determinante \u2248 0): el sistema no tiene una soluci\u00f3n \u00fanica, por eso LU no puede resolverlo. Cambia alg\u00fan coeficiente.',method);return;}}
  else res=mCG(A,b,tol,mx);
  renderSolution(res);
  if(res.type==='lu'){renderLU(res);hideChart();}
  else{renderIterTable(res);showChart(res,method);}
  renderInterp(res,method,A);
  renderPlain(res,method);
}
function syncMethodUI(){
  const m=document.getElementById('m1_method').value;
  document.getElementById('m1_methodDesc').textContent=M1_METHODS[m].desc;
  document.getElementById('m1_omegaWrap').style.display=(m==='sor')?'block':'none';
}
function initM1(){
  buildInputs();syncMethodUI();solve();
  document.getElementById('btnCalc').addEventListener('click',solve);
  document.getElementById('m1_method').addEventListener('change',()=>{syncMethodUI();solve();});
  document.getElementById('m1_omega').addEventListener('input',solve);
  document.getElementById('btnReset').addEventListener('click',()=>{buildInputs();document.getElementById('tol').value='0.01';document.getElementById('maxit').value='100';document.getElementById('m1_method').value='gauss';syncMethodUI();solve();});
  document.getElementById('btnBlock').addEventListener('click',()=>{document.getElementById('a_0_0').value='0.2';solve();});
  document.addEventListener('input',e=>{if(e.target.id&&(e.target.id.startsWith('a_')||e.target.id.startsWith('b_')))solve();});
}

/* arranque */
initM1();
