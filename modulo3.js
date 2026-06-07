/* ===================== MODULO 3 - INTERPOLACION (splines cubicos naturales) ===================== */
const M3_DEF={dias:[1,5,11,16,23,29],pollo:[18,20,25,28,24,21],arroz:[7.0,7.4,8.6,9.2,8.8,8.4],aceite:[14,15,17,18,17,16]};
const M3_PROD=[{key:'pollo',name:'Pollo',unit:'Bs/kg',color:'#E4A64D'},
               {key:'arroz',name:'Arroz',unit:'Bs/kg',color:'#6FC0A8'},
               {key:'aceite',name:'Aceite',unit:'Bs/L',color:'#7FA8D6'}];
let priceChart=null;

function m3BuildGrid(){
  const g=document.getElementById('m3_grid'); g.innerHTML='';
  ['Día','Pollo','Arroz','Aceite'].forEach(h=>{const d=document.createElement('div');d.className='m3-h';d.textContent=h;g.appendChild(d);});
  for(let i=0;i<M3_DEF.dias.length;i++){
    const mk=(id,val,step)=>{const e=document.createElement('input');e.type='number';e.step=step;e.className='num-in';e.id=id;e.value=val;return e;};
    g.appendChild(mk('m3_dia_'+i,M3_DEF.dias[i],'1'));
    g.appendChild(mk('m3_pollo_'+i,M3_DEF.pollo[i],'0.1'));
    g.appendChild(mk('m3_arroz_'+i,M3_DEF.arroz[i],'0.1'));
    g.appendChild(mk('m3_aceite_'+i,M3_DEF.aceite[i],'0.1'));
  }
}
function m3Read(){
  const n=M3_DEF.dias.length,rows=[];
  for(let i=0;i<n;i++) rows.push({x:parseFloat(document.getElementById('m3_dia_'+i).value),
    pollo:parseFloat(document.getElementById('m3_pollo_'+i).value),
    arroz:parseFloat(document.getElementById('m3_arroz_'+i).value),
    aceite:parseFloat(document.getElementById('m3_aceite_'+i).value)});
  rows.sort((a,b)=>a.x-b.x);
  let ok=true; for(let i=1;i<rows.length;i++) if(!(rows[i].x>rows[i-1].x)) ok=false;
  return {rows,ok};
}
// NUCLEO: spline cubico natural (algoritmo de Thomas para el sistema tridiagonal)
function naturalCubicSpline(xs,ys){
  const n=xs.length,a=ys.slice(),b=new Array(n).fill(0),c=new Array(n).fill(0),d=new Array(n).fill(0),h=new Array(n-1);
  for(let i=0;i<n-1;i++) h[i]=xs[i+1]-xs[i];
  const al=new Array(n).fill(0);
  for(let i=1;i<n-1;i++) al[i]=(3/h[i])*(a[i+1]-a[i])-(3/h[i-1])*(a[i]-a[i-1]);
  const l=new Array(n).fill(0),mu=new Array(n).fill(0),z=new Array(n).fill(0);
  l[0]=1;
  for(let i=1;i<n-1;i++){l[i]=2*(xs[i+1]-xs[i-1])-h[i-1]*mu[i-1];mu[i]=h[i]/l[i];z[i]=(al[i]-h[i-1]*z[i-1])/l[i];}
  l[n-1]=1;
  for(let j=n-2;j>=0;j--){c[j]=z[j]-mu[j]*c[j+1];b[j]=(a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;d[j]=(c[j+1]-c[j])/(3*h[j]);}
  return {xs,a,b,c,d,n};
}
function splineEval(sp,x){
  let i=0; const n=sp.n;
  if(x<=sp.xs[0]) i=0; else if(x>=sp.xs[n-1]) i=n-2;
  else { for(let k=0;k<n-1;k++){ if(x>=sp.xs[k]&&x<=sp.xs[k+1]){i=k;break;} } }
  const dx=x-sp.xs[i];
  return sp.a[i]+sp.b[i]*dx+sp.c[i]*dx*dx+sp.d[i]*dx*dx*dx;
}
function m3Splines(data){
  const xs=data.rows.map(r=>r.x),sp={};
  M3_PROD.forEach(p=> sp[p.key]=naturalCubicSpline(xs,data.rows.map(r=>r[p.key])));
  return {sp,xs};
}
function m3RenderChart(data,S){
  const ctx=document.getElementById('priceChart'); if(!ctx) return;
  const x0=S.xs[0],x1=S.xs[S.xs.length-1],ds=[];
  M3_PROD.forEach(p=>{
    const curve=[]; for(let x=x0;x<=x1+1e-9;x+=0.2) curve.push({x:+x.toFixed(2),y:splineEval(S.sp[p.key],x)});
    ds.push({label:p.name,data:curve,borderColor:p.color,backgroundColor:'transparent',borderWidth:2.5,pointRadius:0,tension:0});
    ds.push({label:p.name+' (registros)',data:data.rows.map(r=>({x:r.x,y:r[p.key]})),borderColor:p.color,backgroundColor:p.color,pointRadius:5,pointHoverRadius:6,showLine:false});
  });
  if(priceChart)priceChart.destroy();
  priceChart=new Chart(ctx,{type:'line',data:{datasets:ds},
    options:{responsive:true,parsing:false,interaction:{mode:'nearest',intersect:false},
      plugins:{legend:{labels:{color:'#B6AC94',font:{family:'IBM Plex Mono',size:10},filter:(it)=>!it.text.includes('registros')}}},
      scales:{y:{grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'precio (Bs)',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}},
              x:{type:'linear',grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'día del conflicto',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}}}
    }});
}
function m3RenderEstimate(S){
  const q=Math.max(1,Math.min(30,parseFloat(document.getElementById('m3_query').value)||1));
  document.getElementById('m3_qlabel').textContent=q;
  document.getElementById('m3_estBox').innerHTML=M3_PROD.map(p=>{
    const v=splineEval(S.sp[p.key],q);
    return '<div class="x-row"><div class="name">'+p.name+'<small>estimado el día '+q+'</small></div><div class="qty">Bs '+v.toFixed(2)+' <span>'+p.unit+'</span></div></div>';
  }).join('');
  return q;
}
function m3RenderPlain(S,q){
  const el=document.getElementById('m3_plainConclusion'); if(!el) return;
  const vp=splineEval(S.sp.pollo,q).toFixed(2),va=splineEval(S.sp.arroz,q).toFixed(2),vc=splineEval(S.sp.aceite,q).toFixed(2);
  el.innerHTML='Aunque el día <b>'+q+'</b> nunca se anotó el precio, la curva permite estimarlo: ese día el pollo costaría alrededor de <b>Bs '+vp+'</b>, el arroz <b>Bs '+va+'</b> y el aceite <b>Bs '+vc+'</b>. Así, a partir de solo seis registros, obtenemos una estimación para los 30 días del conflicto.';
}
function m3Calc(){
  const data=m3Read(),el=document.getElementById('m3_plainConclusion');
  if(!data.ok){ document.getElementById('m3_estBox').innerHTML='<p style="color:var(--crisis)" class="mb-0">⚠ Los días deben ser distintos y de menor a mayor. Revisa la columna “Día”.</p>'; if(el) el.textContent=''; return; }
  const S=m3Splines(data); m3RenderChart(data,S); const q=m3RenderEstimate(S); m3RenderPlain(S,q);
}
function initM3(){
  m3BuildGrid(); m3Calc();
  document.getElementById('m3_btnCalc').addEventListener('click',m3Calc);
  document.getElementById('m3_query').addEventListener('input',()=>{const d=m3Read(); if(d.ok){const S=m3Splines(d); const q=m3RenderEstimate(S); m3RenderPlain(S,q);}});
  document.getElementById('m3_btnReset').addEventListener('click',()=>{m3BuildGrid(); document.getElementById('m3_query').value='13'; m3Calc();});
}

/* arranque */
initM3();
