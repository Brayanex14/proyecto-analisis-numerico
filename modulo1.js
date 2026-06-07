/* ===================== MODULO 1 - GAUSS-SEIDEL ===================== */
const DEFAULT_A=[[0.8,0.1,0.1],[0.1,0.7,0.2],[0.1,0.2,0.6]];
const DEFAULT_B=[120,150,180];
const PLANTAS=[{name:"Planta Achacachi",code:"x\u2081 \u00b7 zona norte"},{name:"Planta Senkata",code:"x\u2082 \u00b7 El Alto"},{name:"Planta Patacamaya",code:"x\u2083 \u00b7 zona sur"}];
const HOSP=["H. Municipal Los Andes","Hospital Holand\u00e9s","Hospital La Paz"];
let convChart=null;

function buildInputs(){
  const mg=document.getElementById('matrixGrid'),bg=document.getElementById('bvecGrid');
  mg.innerHTML='';bg.innerHTML='';
  for(let i=0;i<3;i++){
    for(let j=0;j<3;j++){const inp=document.createElement('input');inp.type='number';inp.step='0.01';inp.className='num-in';inp.id=`a_${i}_${j}`;inp.value=DEFAULT_A[i][j];mg.appendChild(inp);}
    const bi=document.createElement('input');bi.type='number';bi.step='1';bi.className='num-in';bi.id=`b_${i}`;bi.value=DEFAULT_B[i];bi.title=HOSP[i];bg.appendChild(bi);
  }
}
function readSystem(){const A=[],b=[];for(let i=0;i<3;i++){A.push([]);for(let j=0;j<3;j++)A[i].push(parseFloat(document.getElementById(`a_${i}_${j}`).value)||0);b.push(parseFloat(document.getElementById(`b_${i}`).value)||0);}return{A,b};}
function checkDominance(A){let dom=true,strict=false;for(let i=0;i<3;i++){const d=Math.abs(A[i][i]);let off=0;for(let j=0;j<3;j++)if(j!==i)off+=Math.abs(A[i][j]);if(d<off)dom=false;if(d>off)strict=true;}return dom&&strict;}
function renderDominance(A){const ok=checkDominance(A);document.getElementById('domBadge').innerHTML=ok?`<span class="dom-badge dom-ok"><span class="dot" style="background:var(--teal)"></span> Diagonalmente dominante \u00b7 convergencia garantizada</span>`:`<span class="dom-badge dom-bad"><span class="dot" style="background:var(--crisis)"></span> NO dominante \u00b7 puede no converger o ser inestable</span>`;return ok;}

// --- NUCLEO: GAUSS-SEIDEL ---
function gaussSeidel(A,b,tolPct,maxIter){
  for(let i=0;i<3;i++)if(A[i][i]===0)return{error:"Hay un cero en la diagonal (a\u1d62\u1d62=0). Reordena las ecuaciones."};
  let x=[0,0,0];const rows=[];let k=0,err=Infinity,diverged=false;
  while(k<maxIter){
    const xOld=[...x];
    for(let i=0;i<3;i++){let s=b[i];for(let j=0;j<3;j++)if(j!==i)s-=A[i][j]*x[j];x[i]=s/A[i][i];}
    err=0;for(let i=0;i<3;i++){const dn=Math.abs(x[i])>1e-12?Math.abs(x[i]):1;err=Math.max(err,Math.abs(x[i]-xOld[i])/dn*100);}
    k++;rows.push({k,x:[...x],err});
    if(!isFinite(x[0])||!isFinite(x[1])||!isFinite(x[2])||err>1e10){diverged=true;break;}
    if(err<tolPct)break;
  }
  return{x,rows,iterations:k,finalErr:err,converged:err<tolPct&&!diverged,diverged};
}
function renderSolution(res){const box=document.getElementById('solBox');if(res.error){box.innerHTML=`<p style="color:var(--crisis)" class="mb-0">\u26a0 ${res.error}</p>`;return;}box.innerHTML=PLANTAS.map((p,i)=>`<div class="x-row"><div class="name">${p.name}<small>${p.code}</small></div><div class="qty">${res.x[i].toFixed(2)} <span>cil.</span></div></div>`).join('');}
function renderIterTable(rows){const b=document.getElementById('iterBody');if(!rows||!rows.length){b.innerHTML=`<tr><td colspan="5">\u2014</td></tr>`;return;}b.innerHTML=rows.map(r=>`<tr><td>${r.k}</td><td>${r.x[0].toFixed(4)}</td><td>${r.x[1].toFixed(4)}</td><td>${r.x[2].toFixed(4)}</td><td>${r.err.toExponential(2)}</td></tr>`).join('');}
function renderChart(rows){const ctx=document.getElementById('convChart');const labels=rows.map(r=>r.k);const data=rows.map(r=>Math.max(r.err,1e-12));if(convChart)convChart.destroy();convChart=new Chart(ctx,{type:'line',data:{labels,datasets:[{label:'Error relativo (%)',data,borderColor:'#E4A64D',backgroundColor:'rgba(228,166,77,.12)',borderWidth:2,pointRadius:3,pointBackgroundColor:'#E4A64D',tension:.25,fill:true}]},options:{responsive:true,plugins:{legend:{labels:{color:'#B6AC94',font:{family:'IBM Plex Mono',size:11}}}},scales:{y:{type:'logarithmic',grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}}},x:{grid:{color:'rgba(255,255,255,.05)'},ticks:{color:'#82785F',font:{family:'IBM Plex Mono',size:10}}}}}});}
function renderInterp(res,dom){const box=document.getElementById('interpBox'),txt=document.getElementById('interpText');if(res.error){box.style.display='none';return;}box.style.display='block';let h='';if(res.converged){h+=`El sistema <b>convergi\u00f3 en ${res.iterations} iteraciones</b> con error final de <b>${res.finalErr.toExponential(2)} %</b>. La distribuci\u00f3n asigna <b>${res.x[0].toFixed(1)}</b>, <b>${res.x[1].toFixed(1)}</b> y <b>${res.x[2].toFixed(1)}</b> cilindros a Achacachi, Senkata y Patacamaya. `;h+=dom?`Como la matriz es diagonalmente dominante, la convergencia estaba garantizada y fue estable.`:`Pese a no ser estrictamente dominante, esta configuraci\u00f3n a\u00fan converge \u2014 la dominancia es <i>suficiente</i> pero no <i>necesaria</i>.`;}else if(res.diverged){h+=`<b style="color:var(--crisis)">El sistema divergi\u00f3.</b> Al perder la dominancia diagonal (p. ej. un bloqueo que reduce una ruta principal), Gauss-Seidel ya no se estabiliza: las correcciones crecen sin control. Es lo que ocurre cuando un corte de ruta desarticula la red.`;}else{h+=`No se alcanz\u00f3 la tolerancia en ${res.iterations} iteraciones (error <b>${res.finalErr.toExponential(2)} %</b>). El sistema converge muy lento: sube el m\u00e1ximo de iteraciones o revisa la dominancia.`;}txt.innerHTML=h;}
function renderPlain(res){
  const el=document.getElementById('plainConclusion'); if(!el) return;
  if(res.error){ el.innerHTML='Hay un dato que impide resolver el sistema. Revisa los valores ingresados.'; return; }
  if(res.converged){
    el.innerHTML=`Con los datos actuales, para que los tres hospitales reciban el oxígeno que necesitan, las plantas deberían despachar aproximadamente <b>${Math.round(res.x[0])} cilindros</b> desde Achacachi, <b>${Math.round(res.x[1])}</b> desde Senkata y <b>${Math.round(res.x[2])}</b> desde Patacamaya cada día. La computadora encontró este reparto tras <b>${res.iterations} ajustes</b>, cuando los números dejaron de cambiar.`;
  } else if(res.diverged){
    el.innerHTML=`Con esta configuración el método <b>no logra encontrar una respuesta</b>: en vez de estabilizarse, los números se disparan. Significa que, tal como quedó planteada la red de rutas, <b>el oxígeno no se puede repartir de forma que alcance para todos los hospitales</b>.`;
  } else {
    el.innerHTML=`El método se quedó sin intentos antes de afinar la respuesta (alcanzó el máximo de ajustes). Todavía no es una solución confiable: habría que darle más iteraciones o revisar la red de rutas.`;
  }
}
function calcular(){const{A,b}=readSystem();const dom=renderDominance(A);const tol=parseFloat(document.getElementById('tol').value)||0.01;const mx=parseInt(document.getElementById('maxit').value)||100;const res=gaussSeidel(A,b,tol,mx);renderSolution(res);renderIterTable(res.rows);if(res.rows)renderChart(res.rows);renderInterp(res,dom);renderPlain(res);}
function initM1(){
  buildInputs();renderDominance(DEFAULT_A);calcular();
  document.getElementById('btnCalc').addEventListener('click',calcular);
  document.getElementById('btnReset').addEventListener('click',()=>{buildInputs();document.getElementById('tol').value='0.01';document.getElementById('maxit').value='100';calcular();});
  document.getElementById('btnBlock').addEventListener('click',()=>{document.getElementById('a_0_0').value='0.2';renderDominance(readSystem().A);calcular();});
  document.addEventListener('input',e=>{if(e.target.id&&e.target.id.startsWith('a_'))renderDominance(readSystem().A);});
}

/* arranque */
initM1();
