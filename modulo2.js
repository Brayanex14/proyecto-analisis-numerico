/* ===================== MODULO 2 - EDO RK4 (vaciado de silo) ===================== */
let siloChart=null; let m2CompareChart=null;
function m2val(id){return parseFloat(document.getElementById(id).value);}
function m2Read(){
  return {R0:m2val('m2_R0')||0, C:m2val('m2_C')||0, E:m2val('m2_E')||0,
          h:Math.max(m2val('m2_h')||0.5,0.001), tmax:Math.max(m2val('m2_tmax')||15,0.1)};
}
// NUCLEO: Runge-Kutta 4 sobre R'(t)=E-C
function rk4Silo(p){
  const f=(t,R)=> p.E - p.C;
  let t=0,R=p.R0,n=0; const rows=[]; const series=[{x:0,y:R}]; let emptyDay=null;
  const net=p.C-p.E;                       // tasa neta de salida (t/dia)
  while(t < p.tmax-1e-9 && n < 100000){
    const k1=f(t,R), k2=f(t+p.h/2,R+p.h/2*k1), k3=f(t+p.h/2,R+p.h/2*k2), k4=f(t+p.h,R+p.h*k3);
    let Rn=R + p.h/6*(k1+2*k2+2*k3+k4);
    let tn=t+p.h;
    if(Rn<=0 && net>0){                    // cruce exacto por cero
      const tc=t + R/net;
      rows.push({n:n+1,t:tc,k1,k2,k3,k4,R:0});
      series.push({x:tc,y:0}); emptyDay=tc; break;
    }
    rows.push({n:n+1,t:tn,k1,k2,k3,k4,R:Rn});
    series.push({x:tn,y:Rn});
    R=Rn; t=tn; n++;
  }
  return {rows,series,emptyDay,net,p};
}
function m2RenderHead(res){
  const box=document.getElementById('m2_emptyBox'), p=res.p;
  if(res.net<=0){
    box.innerHTML=`<div class="x-row"><div class="name">El silo NO se vac\u00eda<small>la entrada iguala o supera al consumo</small></div><div class="qty" style="color:var(--teal)">\u221e</div></div>`;
    return;
  }
  const d=res.emptyDay!=null?res.emptyDay:(p.R0/res.net);
  const dias=Math.floor(d), hrs=Math.round((d-dias)*24);
  box.innerHTML=`<div class="x-row"><div class="name">Autonom\u00eda del silo<small>R(0)=${p.R0} t \u00b7 salida neta ${res.net.toFixed(2)} t/d\u00eda</small></div><div class="qty">${d.toFixed(2)} <span>d\u00edas</span></div></div>`+
  `<div class="x-row"><div class="name">Equivale a<small>tiempo hasta reserva cero</small></div><div class="qty" style="font-size:1.15rem;color:var(--text-dim)">${dias} d ${hrs} h</div></div>`;
}
function m2RenderTable(rows){
  document.getElementById('m2_iterBody').innerHTML=rows.map(r=>`<tr><td>${r.n}</td><td>${r.t.toFixed(3)}</td><td>${r.k1.toFixed(3)}</td><td>${r.k2.toFixed(3)}</td><td>${r.k3.toFixed(3)}</td><td>${r.k4.toFixed(3)}</td><td>${r.R.toFixed(4)}</td></tr>`).join('');
}
function m2RenderChart(res){
  const ctx=document.getElementById('siloChart');
  if(siloChart)siloChart.destroy();
  siloChart=new Chart(ctx,{type:'line',
    data:{datasets:[{label:'Reserva R(t) \u2014 t',data:res.series,borderColor:'#E4A64D',backgroundColor:'rgba(228,166,77,.14)',borderWidth:2.5,pointRadius:2,pointBackgroundColor:'#E4A64D',tension:0,fill:true}]},
    options:{responsive:true,parsing:false,
      plugins:{legend:{labels:{color:'#B6AC94',font:{family:'IBM Plex Mono',size:11}}}},
      scales:{y:{min:0,grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'reserva (t)',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}},
              x:{type:'linear',grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'tiempo (d\u00edas)',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}}}
    }});
}
function m2RenderInterp(res){
  const box=document.getElementById('m2_interpBox'), txt=document.getElementById('m2_interpText'), p=res.p;
  box.style.display='block'; let h='';
  if(res.net<=0){
    h=`Con un consumo de <b>${p.C} t/d\u00eda</b> y una entrada de <b>${p.E} t/d\u00eda</b>, la tasa neta es \u2265 0: la reserva se mantiene o crece y <b>el silo nunca se vac\u00eda</b>. En el conflicto, significa que el abastecimiento alterno alcanza para cubrir el consumo \u2014 el bloqueo deja de ser cr\u00edtico.`;
  } else {
    const an=p.R0/res.net;
    h=`Con la reserva inicial de <b>${p.R0} t</b> y una salida neta de <b>${res.net.toFixed(2)} t/d\u00eda</b> (consumo ${p.C} \u2212 entrada ${p.E}), el silo se agota en <b>${res.emptyDay.toFixed(2)} d\u00edas</b>. Como R\u2032 es constante, la soluci\u00f3n exacta es la recta R(t)=${p.R0}\u2212${res.net.toFixed(2)}\u00b7t, que predice el vaciado en <b>${an.toFixed(2)} d\u00edas</b>: RK4 reproduce el valor anal\u00edtico sin error apreciable, lo que <b>valida la implementaci\u00f3n</b>. Cada d\u00eda de bloqueo acerca a la granja al punto en que las 30.000 aves se quedan sin alimento.`;
  }
  txt.innerHTML=h;
}
function m2RenderCompare(p){
  const ctx=document.getElementById('siloCompareChart'); if(!ctx) return;
  const C=p.C, R0=p.R0;
  const tBase = C>0 ? R0/C : 15;
  const tComp = Math.max(p.tmax, tBase*1.5, 1);
  const N=60;
  function lineFor(E){const net=C-E;const a=[];for(let i=0;i<=N;i++){const t=tComp*i/N;a.push({x:+t.toFixed(3),y:Math.max(R0-net*t,0)});}return a;}
  const eMid=+(C/2).toFixed(2);
  const ds=[
    {label:'Bloqueo total (E=0)',data:lineFor(0),borderColor:'#E05C3C',backgroundColor:'transparent',borderWidth:2.5,pointRadius:0,tension:0},
    {label:'Goteo parcial (E='+eMid+')',data:lineFor(eMid),borderColor:'#E4A64D',backgroundColor:'transparent',borderWidth:2.5,pointRadius:0,tension:0},
    {label:'Entrada que empata (E='+C+')',data:lineFor(C),borderColor:'#6FC0A8',backgroundColor:'transparent',borderWidth:2.5,pointRadius:0,borderDash:[6,5],tension:0}
  ];
  if(m2CompareChart)m2CompareChart.destroy();
  m2CompareChart=new Chart(ctx,{type:'line',data:{datasets:ds},
    options:{responsive:true,parsing:false,
      plugins:{legend:{labels:{color:'#B6AC94',font:{family:'IBM Plex Mono',size:10}}}},
      scales:{y:{min:0,grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'reserva (t)',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}},
              x:{type:'linear',grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}},title:{display:true,text:'tiempo (d\u00edas)',color:'#82785F',font:{family:'IBM Plex Mono',size:10}}}}
    }});
}
function m2RenderPlain(res){
  const el=document.getElementById('m2_plainConclusion'); if(!el) return; const p=res.p;
  if(res.net<=0){
    el.innerHTML='Con la entrada actual (<b>'+p.E+' t/d\u00eda</b>), el alimento que llega alcanza para cubrir lo que comen las aves: <b>la reserva no se agota</b>. Para esta granja, el bloqueo deja de ser un problema de alimento.';
  } else {
    const d=res.emptyDay!=null?res.emptyDay:(p.R0/res.net); const dias=Math.floor(d);
    el.innerHTML='Con '+p.R0+' toneladas guardadas y un consumo neto de <b>'+res.net.toFixed(2)+' t/d\u00eda</b>, el silo alcanza para unos <b>'+d.toFixed(1)+' d\u00edas</b> (cerca de '+dias+' d\u00edas). Pasado ese punto, las 30.000 aves se quedan sin alimento. Por eso lograr que entre aunque sea <b>un poco</b> de comida por rutas alternas estira mucho ese plazo \u2014como se ve en la gr\u00e1fica de comparaci\u00f3n\u2014.';
  }
}
function m2Calc(){const p=m2Read();const res=rk4Silo(p);m2RenderHead(res);m2RenderTable(res.rows);m2RenderChart(res);m2RenderInterp(res);m2RenderCompare(res.p);m2RenderPlain(res);}
function initM2(){
  m2Calc();
  document.getElementById('m2_btnCalc').addEventListener('click',m2Calc);
  document.getElementById('m2_btnDrip').addEventListener('click',()=>{document.getElementById('m2_E').value='0.5';m2Calc();});
  document.getElementById('m2_btnReset').addEventListener('click',()=>{document.getElementById('m2_R0').value='15';document.getElementById('m2_C').value='2.1';document.getElementById('m2_E').value='0';document.getElementById('m2_h').value='0.5';document.getElementById('m2_tmax').value='15';m2Calc();});
}

/* arranque */
initM2();
