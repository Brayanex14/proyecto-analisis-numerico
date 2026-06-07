/* ===================== MODULO 5 - RAICES DE ECUACIONES (Newton-Raphson, Biseccion, Secante) ===================== */
const M5_DEF={dias:[1,5,11,16,23,29],pollo:[18,20,25,28,24,21],arroz:[7.0,7.4,8.6,9.2,8.8,8.4],aceite:[14,15,17,18,17,16]};
const M5_Q={pollo:1,arroz:0.6,aceite:0.15};
let rootChart=null;

// --- spline cubico natural (mismo nucleo) ---
function naturalCubicSpline(xs,ys){
  const n=xs.length,a=ys.slice(),b=new Array(n).fill(0),c=new Array(n).fill(0),d=new Array(n).fill(0),h=new Array(n-1);
  for(let i=0;i<n-1;i++) h[i]=xs[i+1]-xs[i];
  const al=new Array(n).fill(0);
  for(let i=1;i<n-1;i++) al[i]=(3/h[i])*(a[i+1]-a[i])-(3/h[i-1])*(a[i]-a[i-1]);
  const l=new Array(n).fill(0),mu=new Array(n).fill(0),z=new Array(n).fill(0); l[0]=1;
  for(let i=1;i<n-1;i++){l[i]=2*(xs[i+1]-xs[i-1])-h[i-1]*mu[i-1];mu[i]=h[i]/l[i];z[i]=(al[i]-h[i-1]*z[i-1])/l[i];}
  l[n-1]=1;
  for(let j=n-2;j>=0;j--){c[j]=z[j]-mu[j]*c[j+1];b[j]=(a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;d[j]=(c[j+1]-c[j])/(3*h[j]);}
  return {xs,a,b,c,d,n};
}
function splineEval(sp,x){
  let i=0;const n=sp.n;
  if(x<=sp.xs[0])i=0; else if(x>=sp.xs[n-1])i=n-2;
  else { for(let k=0;k<n-1;k++){ if(x>=sp.xs[k]&&x<=sp.xs[k+1]){i=k;break;} } }
  const dx=x-sp.xs[i]; return sp.a[i]+sp.b[i]*dx+sp.c[i]*dx*dx+sp.d[i]*dx*dx*dx;
}
// integral de un spline de t0 a t (Simpson fino)
function integTo(sp,t0,t,n){ if(t<=t0) return 0; n=Math.max(4,n); if(n%2)n++; const h=(t-t0)/n; let s=splineEval(sp,t0)+splineEval(sp,t); for(let i=1;i<n;i++) s+=(i%2?4:2)*splineEval(sp,t0+i*h); return s*h/3; }

// --- construir f(t) y f'(t) ---
function buildFun(B){
  const xs=M5_DEF.dias, t0=xs[0], t1=xs[xs.length-1];
  const sp={pollo:naturalCubicSpline(xs,M5_DEF.pollo),arroz:naturalCubicSpline(xs,M5_DEF.arroz),aceite:naturalCubicSpline(xs,M5_DEF.aceite)};
  const G=(t)=> M5_Q.pollo*integTo(sp.pollo,t0,t,80)+M5_Q.arroz*integTo(sp.arroz,t0,t,80)+M5_Q.aceite*integTo(sp.aceite,t0,t,80);
  const f=(t)=> G(t)-B;
  const fp=(t)=> M5_Q.pollo*splineEval(sp.pollo,t)+M5_Q.arroz*splineEval(sp.arroz,t)+M5_Q.aceite*splineEval(sp.aceite,t);
  return {f,fp,G,t0,t1};
}

// --- METODOS ---
function newton(f,fp,x0,tol,maxit){
  const rows=[]; let x=x0;
  for(let k=0;k<maxit;k++){ const fx=f(x),d=fp(x); if(Math.abs(d)<1e-14)return{rows,diverged:true}; const xn=x-fx/d; rows.push({k:k+1,x:xn,fx:f(xn),err:Math.abs(xn-x)}); if(Math.abs(xn-x)<tol)return{root:xn,rows,iterations:k+1,converged:true}; x=xn; }
  return{root:x,rows,iterations:maxit,converged:false};
}
function biseccion(f,a,b,tol,maxit){
  let fa=f(a),fb=f(b); if(fa*fb>0)return{noRoot:true}; const rows=[]; let c=a;
  for(let k=0;k<maxit;k++){ c=(a+b)/2; const fc=f(c); rows.push({k:k+1,x:c,fx:fc,err:(b-a)/2}); if(Math.abs(fc)<1e-12||(b-a)/2<tol)return{root:c,rows,iterations:k+1,converged:true}; if(fa*fc<0){b=c;fb=fc;}else{a=c;fa=fc;} }
  return{root:c,rows,iterations:maxit,converged:false};
}
function secante(f,x0,x1,tol,maxit){
  const rows=[]; let f0=f(x0),f1=f(x1);
  for(let k=0;k<maxit;k++){ if(Math.abs(f1-f0)<1e-14)return{rows,diverged:true}; const x2=x1-f1*(x1-x0)/(f1-f0); rows.push({k:k+1,x:x2,fx:f(x2),err:Math.abs(x2-x1)}); if(Math.abs(x2-x1)<tol)return{root:x2,rows,iterations:k+1,converged:true}; x0=x1;f0=f1;x1=x2;f1=f(x2); }
  return{root:x1,rows,iterations:maxit,converged:false};
}
function estOrder(rows,r){
  const e=rows.map(row=>Math.abs(row.x-r)).filter(v=>v>1e-13);
  if(e.length<3)return null; const n=e.length;
  const p=Math.log(e[n-1]/e[n-2])/Math.log(e[n-2]/e[n-3]);
  return (isFinite(p)&&p>0&&p<5)?p:null;
}

// --- render ---
function m5RenderRes(root,B,G,t1){
  const total=G(t1); const box=document.getElementById('m5_resBox');
  const com=v=>v.toFixed(2).replace('.',',');
  const dia=Math.round(root), hrs=Math.round((root-Math.floor(root))*24);
  box.innerHTML=`<div class="x-row"><div class="name">Día cero<small>cuando el gasto acumulado = Bs ${B}</small></div><div class="qty">\u2248 día ${dia}</div></div>`+
   `<div class="x-row"><div class="name">Valor exacto del método<small>t = ${com(root)} días \u00b7 el decimal es una fracción del día (\u2248 ${hrs} h), no una fecha</small></div><div class="qty" style="font-size:1.05rem;color:var(--text-dim)">${com(root)}</div></div>`+
   `<div class="x-row"><div class="name">Gasto total del mes<small>si comprara los 29 días</small></div><div class="qty" style="font-size:1.15rem;color:var(--teal)">Bs ${total.toFixed(0)}</div></div>`;
}
function m5RenderChart(F,root){
  const ctx=document.getElementById('rootChart'); if(!ctx)return;
  const curve=[]; for(let t=F.t0;t<=F.t1+1e-9;t+=0.2) curve.push({x:+t.toFixed(2),y:F.f(t)});
  const ds=[
    {label:'f(t) = gasto acumulado − presupuesto',data:curve,borderColor:'#E4A64D',backgroundColor:'rgba(228,166,77,.10)',borderWidth:2.5,pointRadius:0,fill:'origin',tension:0},
    {label:'cero',data:[{x:F.t0,y:0},{x:F.t1,y:0}],borderColor:'#9aa0aa',borderDash:[6,5],borderWidth:1.5,pointRadius:0,fill:false}
  ];
  if(root!=null) ds.push({label:'día cero',data:[{x:root,y:0}],borderColor:'#E05C3C',backgroundColor:'#E05C3C',pointRadius:7,pointHoverRadius:8,showLine:false});
  if(rootChart)rootChart.destroy();
  rootChart=new Chart(ctx,{type:'line',data:{datasets:ds},
    options:{responsive:true,parsing:false,
      plugins:{legend:{labels:{color:'#B6AC94',font:{family:'IBM Plex Mono',size:10}}}},
      scales:{y:{grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'f(t) en Bs',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}},
              x:{type:'linear',grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'día del conflicto',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}}}
    }});
}
function m5RenderTable(run,method){
  const names={newton:'Newton-Raphson',biseccion:'Bisección',secante:'Secante'};
  document.getElementById('m5_detailLabel').textContent='Tabla de iteraciones · '+names[method];
  if(!run.rows||!run.rows.length){document.getElementById('m5_detail').innerHTML='<p class="dim" style="padding:.8rem">Sin iteraciones.</p>';return;}
  let h='<table class="tbl-iter"><thead><tr><th>k</th><th>t (día)</th><th>f(t)</th><th>error</th></tr></thead><tbody>';
  h+=run.rows.map(r=>`<tr><td>${r.k}</td><td>${r.x.toFixed(6)}</td><td>${r.fx.toFixed(4)}</td><td>${r.err.toExponential(2)}</td></tr>`).join('');
  document.getElementById('m5_detail').innerHTML=h+'</tbody></table>';
}
function m5RenderCmp(runs,r){
  const names={newton:'Newton-Raphson',biseccion:'Bisección',secante:'Secante'};
  const teor={newton:'2',biseccion:'1',secante:'1.6'};
  document.getElementById('m5_cmpBody').innerHTML=['newton','secante','biseccion'].map(k=>{
    const run=runs[k]; const it=run.iterations!=null?run.iterations:'—';
    const p=run.rows?estOrder(run.rows,r):null;
    const pTxt=p!=null?p.toFixed(2):'—';
    return `<tr><td>${names[k]}</td><td>${it}</td><td>${pTxt} <span style="color:var(--text-faint)">(teórico ≈ ${teor[k]})</span></td></tr>`;
  }).join('');
}
function m5RenderInterp(runs,method,r,F,B){
  const box=document.getElementById('m5_interpBox'),txt=document.getElementById('m5_interpText'); box.style.display='block';
  const it=runs[method].iterations, names={newton:'Newton-Raphson',biseccion:'Bisección',secante:'Secante'};
  txt.innerHTML=`Con <b>${names[method]}</b> se encontró que el gasto acumulado alcanza el presupuesto de <b>Bs ${B}</b> alrededor del <b>día ${Math.round(r)}</b> (t = ${r.toFixed(2).replace('.',',')} días), en <b>${it}</b> iteraciones. Comparando los tres métodos: <b>Newton</b> es el más rápido (orden ≈ 2) pero exige una buena estimación inicial; la <b>Secante</b> es casi tan rápida sin usar la derivada; y la <b>Bisección</b> es la más lenta pero la más robusta, porque solo necesita que la función cambie de signo en el intervalo.`;
}
function m5RenderPlain(r,B){
  const el=document.getElementById('m5_plainConclusion'); if(!el)return;
  const dia=Math.round(r); const com=r.toFixed(2).replace('.',',');
  el.innerHTML=`Con un presupuesto de <b>Bs ${B}</b> para pollo, arroz y aceite, la familia agotaría ese dinero <b>alrededor del día ${dia}</b> del mes (el método lo ubica con precisión en t = ${com} días, donde el decimal es una fracción del día, no una fecha). A partir de ahí ya no le alcanzaría para mantener el mismo consumo: tendría que comprar menos o buscar productos más baratos.`;
}
function m5NoRoot(F,B,total){
  document.getElementById('m5_resBox').innerHTML=`<div class="blocked-msg">⚠ No hay “día cero” este mes<p>Con un presupuesto de Bs ${B}, el gasto del mes (≈ Bs ${total.toFixed(0)}) ${B<=0?'ya estaría superado desde el inicio':'nunca llega a ese límite'}. ${B>0?'Baja el presupuesto por debajo de Bs '+total.toFixed(0)+' para ver un día cero.':''}</p></div>`;
  m5RenderChart(F,null);
  document.getElementById('m5_detail').innerHTML='<p class="dim" style="padding:.8rem">No se puede iterar: la función no cruza el cero dentro del mes.</p>';
  document.getElementById('m5_cmpBody').innerHTML='<tr><td colspan="3" style="text-align:center;color:var(--text-faint)">—</td></tr>';
  const box=document.getElementById('m5_interpBox');box.style.display='block';
  document.getElementById('m5_interpText').innerHTML=`La función f(t)=gasto−presupuesto no cambia de signo en el mes, así que no tiene raíz en ese intervalo. Es un resultado válido: significa que el presupuesto alcanza para todo el mes (o que ya estaba agotado de entrada).`;
  document.getElementById('m5_plainConclusion').innerHTML=`Con ese presupuesto la familia <b>no se queda sin fondos</b> dentro del mes para estos tres productos. Prueba con un presupuesto menor para encontrar el día crítico.`;
}
function m5Calc(){
  const B=parseFloat(document.getElementById('m5_budget').value)||0;
  const x0=parseFloat(document.getElementById('m5_x0').value)||15;
  const tol=parseFloat(document.getElementById('m5_tol').value)||1e-4;
  const maxit=parseInt(document.getElementById('m5_maxit').value)||50;
  const method=document.getElementById('m5_method').value;
  const F=buildFun(B); const total=F.G(F.t1);
  if(B<=0 || B>=total){ m5NoRoot(F,B,total); return; }
  const ref=biseccion(F.f,F.t0,F.t1,1e-11,300); const r=ref.root;
  const runs={
    newton:newton(F.f,F.fp,x0,tol,maxit),
    biseccion:biseccion(F.f,F.t0,F.t1,tol,maxit),
    secante:secante(F.f,x0,Math.min(x0+4,F.t1),tol,maxit)
  };
  const sel=runs[method]; const root=(sel&&sel.root!=null)?sel.root:r;
  m5RenderRes(root,B,F.G,F.t1);
  m5RenderChart(F,root);
  m5RenderTable(sel,method);
  m5RenderCmp(runs,r);
  m5RenderInterp(runs,method,r,F,B);
  m5RenderPlain(root,B);
}
function initM5(){
  m5Calc();
  document.getElementById('m5_btnCalc').addEventListener('click',m5Calc);
  document.getElementById('m5_method').addEventListener('change',m5Calc);
  document.getElementById('m5_btnReset').addEventListener('click',()=>{document.getElementById('m5_budget').value='600';document.getElementById('m5_x0').value='15';document.getElementById('m5_tol').value='0.0001';document.getElementById('m5_maxit').value='50';document.getElementById('m5_method').value='newton';m5Calc();});
}

/* arranque */
initM5();
