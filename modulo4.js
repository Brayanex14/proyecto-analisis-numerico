/* ===================== MODULO 4 - INTEGRACION NUMERICA (Trapecio, Simpson 1/3, Simpson 3/8) ===================== */
const M4_DEF={dias:[1,5,11,16,23,29],pollo:[18,20,25,28,24,21],arroz:[7.0,7.4,8.6,9.2,8.8,8.4],aceite:[14,15,17,18,17,16]};
const M4_PROD=[{key:'pollo',name:'Pollo',unit:'kg',qid:'m4_q_pollo',color:'#E4A64D',fill:'rgba(228,166,77,.18)'},{key:'arroz',name:'Arroz',unit:'kg',qid:'m4_q_arroz',color:'#6FC0A8',fill:'rgba(111,192,168,.18)'},{key:'aceite',name:'Aceite',unit:'L',qid:'m4_q_aceite',color:'#7FA8D6',fill:'rgba(127,168,214,.18)'}];
let costChart=null; let areaChart=null; let m4Last=null;

// --- spline cubico natural (mismo nucleo del Modulo 3) ---
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
// --- muestreo y reglas de integracion ---
function sample(sp,t0,t1,n){const h=(t1-t0)/n,f=[];for(let i=0;i<=n;i++)f.push(splineEval(sp,t0+i*h));return {f,h};}
function trapecio(f,h){let s=f[0]+f[f.length-1];for(let i=1;i<f.length-1;i++)s+=2*f[i];return s*h/2;}
function simpson13(f,h){const n=f.length-1;let s=f[0]+f[n];for(let i=1;i<n;i++)s+=(i%2?4:2)*f[i];return s*h/3;}
function simpson38(f,h){const n=f.length-1;let s=f[0]+f[n];for(let i=1;i<n;i++)s+=(i%3===0?2:3)*f[i];return s*3*h/8;}

function m4ReadN(){let n=Math.round((parseFloat(document.getElementById('m4_n').value)||30)/6)*6; if(n<6)n=6; document.getElementById('m4_n').value=n; return n;}
function m4Calc(){
  const n=m4ReadN();
  const xs=M4_DEF.dias, t0=xs[0], t1=xs[xs.length-1], span=t1-t0;
  let per=[], totTrap=0,totS13=0,totS38=0,totRef=0,totNormal=0;
  M4_PROD.forEach(p=>{
    const q=parseFloat(document.getElementById(p.qid).value)||0;
    const sp=naturalCubicSpline(xs,M4_DEF[p.key]);
    const s=sample(sp,t0,t1,n), sref=sample(sp,t0,t1,600);
    const iTrap=trapecio(s.f,s.h), iS13=simpson13(s.f,s.h), iS38=simpson38(s.f,s.h), iRef=simpson13(sref.f,sref.h);
    const gNormal=q*M4_DEF[p.key][0]*span;   // precio del dia 1 mantenido
    const item={key:p.key,name:p.name,color:p.color,q,sp,integral:iS13,normalPrice:M4_DEF[p.key][0],gS13:q*iS13,gRef:q*iRef,gNormal,perdida:q*iRef-gNormal};
    per.push(item);
    totTrap+=q*iTrap; totS13+=q*iS13; totS38+=q*iS38; totRef+=q*iRef; totNormal+=gNormal;
  });
  const perdida=totRef-totNormal;
  m4Last={t0,t1,items:{}}; per.forEach(it=>m4Last.items[it.key]=it);
  m4RenderArea();
  m4RenderRes(totS13,totNormal,perdida);
  m4RenderChart(per);
  m4RenderMethods({totTrap,totS13,totS38,totRef});
  m4RenderInterp(per,totRef,totNormal,perdida);
  m4RenderPlain(per,totRef,totNormal,perdida);
}
function m4RenderRes(crisis,normal,perdida){
  document.getElementById('m4_resBox').innerHTML=
   `<div class="x-row"><div class="name">Gasto del mes (en crisis)<small>pollo + arroz + aceite</small></div><div class="qty">Bs ${crisis.toFixed(0)}</div></div>`+
   `<div class="x-row"><div class="name">Si los precios no sub\u00edan<small>al precio del d\u00eda 1</small></div><div class="qty" style="color:var(--teal)">Bs ${normal.toFixed(0)}</div></div>`+
   `<div class="x-row"><div class="name">P\u00e9rdida del poder adquisitivo<small>lo que se pag\u00f3 de m\u00e1s</small></div><div class="qty" style="color:var(--crisis)">Bs ${perdida.toFixed(0)}</div></div>`;
}
function m4RenderChart(per){
  const ctx=document.getElementById('costChart'); if(!ctx) return;
  if(costChart)costChart.destroy();
  costChart=new Chart(ctx,{type:'bar',
    data:{labels:per.map(p=>p.name),datasets:[
      {label:'En crisis',data:per.map(p=>+p.gRef.toFixed(1)),backgroundColor:'#E05C3C'},
      {label:'Si no sub\u00edan',data:per.map(p=>+p.gNormal.toFixed(1)),backgroundColor:'#6FC0A8'}
    ]},
    options:{responsive:true,
      plugins:{legend:{labels:{color:'#B6AC94',font:{family:'IBM Plex Mono',size:10}}}},
      scales:{y:{grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'gasto del mes (Bs)',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}},
              x:{grid:{display:false},ticks:{color:'#B6AC94',font:{family:'IBM Plex Mono',size:11}}}}
    }});
}
function m4RenderMethods(m){
  const rows=[['Trapecio',m.totTrap],['Simpson 1/3',m.totS13],['Simpson 3/8',m.totS38]];
  document.getElementById('m4_methodBody').innerHTML=rows.map(r=>{
    const err=Math.abs(r[1]-m.totRef);
    return `<tr><td>${r[0]}</td><td>${r[1].toFixed(3)}</td><td>${err.toExponential(2)}</td></tr>`;
  }).join('')+`<tr><td>Referencia (n=600)</td><td>${m.totRef.toFixed(3)}</td><td>\u2014</td></tr>`;
  const errT=Math.abs(m.totTrap-m.totRef),e13=Math.abs(m.totS13-m.totRef),e38=Math.abs(m.totS38-m.totRef);
  const best=Math.min(errT,e13,e38);
  const name=best===e13?'Simpson 1/3':(best===e38?'Simpson 3/8':'Trapecio');
  document.getElementById('m4_methodNote').innerHTML=`El m\u00e9todo m\u00e1s preciso (menor error frente a la referencia fina) fue <b style="color:var(--text)">${name}</b>. Los m\u00e9todos de Simpson ganan al Trapecio porque ajustan la curva con par\u00e1bolas.`;
}
function m4RenderInterp(per,crisis,normal,perdida){
  const box=document.getElementById('m4_interpBox'),txt=document.getElementById('m4_interpText');
  box.style.display='block';
  const top=per.slice().sort((a,b)=>b.perdida-a.perdida)[0];
  const pct=normal>0?(perdida/normal*100):0;
  txt.innerHTML=`Durante el mes la familia habr\u00eda gastado <b>Bs ${crisis.toFixed(0)}</b> en estos tres productos. A precios del d\u00eda 1 habr\u00eda gastado <b>Bs ${normal.toFixed(0)}</b>; la subida le cost\u00f3 <b>Bs ${perdida.toFixed(0)}</b> extra, una p\u00e9rdida del <b>${pct.toFixed(1)}%</b> de poder adquisitivo. El producto que m\u00e1s pes\u00f3 en ese aumento fue el <b>${top.name}</b>.`;
}
function m4RenderPlain(per,crisis,normal,perdida){
  const el=document.getElementById('m4_plainConclusion'); if(!el) return;
  const top=per.slice().sort((a,b)=>b.perdida-a.perdida)[0];
  el.innerHTML=`Sumando el costo de cada d\u00eda, la familia gast\u00f3 alrededor de <b>Bs ${crisis.toFixed(0)}</b> en pollo, arroz y aceite. Si los precios se hubieran mantenido como al inicio, habr\u00eda gastado <b>Bs ${normal.toFixed(0)}</b>. Esos <b>Bs ${perdida.toFixed(0)}</b> de diferencia son el golpe real del conflicto al bolsillo, y el <b>${top.name}</b> fue el que m\u00e1s lo encareci\u00f3.`;
}
function m4RenderArea(){
  if(!m4Last) return;
  const sel=document.getElementById('m4_prodSel'); if(!sel) return;
  const key=sel.value, prod=M4_PROD.find(p=>p.key===key), it=m4Last.items[key], t0=m4Last.t0, t1=m4Last.t1;
  const curve=[]; for(let x=t0;x<=t1+1e-9;x+=0.2) curve.push({x:+x.toFixed(2),y:splineEval(it.sp,x)});
  const base=[{x:t0,y:it.normalPrice},{x:t1,y:it.normalPrice}];
  if(areaChart)areaChart.destroy();
  areaChart=new Chart(document.getElementById('areaChart'),{type:'line',
    data:{datasets:[
      {label:'Precio '+prod.name+' P(t)',data:curve,borderColor:prod.color,backgroundColor:prod.fill,borderWidth:2.5,pointRadius:0,fill:'origin',tension:0},
      {label:'Precio del d\u00eda 1',data:base,borderColor:'#9aa0aa',borderDash:[6,5],borderWidth:1.5,pointRadius:0,fill:false}
    ]},
    options:{responsive:true,parsing:false,
      plugins:{legend:{labels:{color:'#B6AC94',font:{family:'IBM Plex Mono',size:10}}}},
      scales:{y:{min:0,grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'precio (Bs/'+prod.unit+')',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}},
              x:{type:'linear',grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'d\u00eda del conflicto',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}}}
    }});
  document.getElementById('m4_areaNote').innerHTML='El \u00e1rea sombreada bajo la curva es la integral del precio del '+prod.name+' (\u2248 '+it.integral.toFixed(1)+' Bs\u00b7d\u00eda). Al multiplicarla por el consumo diario ('+it.q+' '+prod.unit+'/d\u00eda) se obtiene el gasto del mes en '+prod.name+': <b>Bs '+(it.integral*it.q).toFixed(0)+'</b>.';
}
function initM4(){
  m4Calc();
  document.getElementById('m4_prodSel').addEventListener('change',m4RenderArea);
  document.getElementById('m4_btnCalc').addEventListener('click',m4Calc);
  document.getElementById('m4_btnReset').addEventListener('click',()=>{document.getElementById('m4_q_pollo').value='1';document.getElementById('m4_q_arroz').value='0.6';document.getElementById('m4_q_aceite').value='0.15';document.getElementById('m4_n').value='30';m4Calc();});
}

/* arranque */
initM4();
