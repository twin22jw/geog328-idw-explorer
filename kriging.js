(() => {
'use strict';
const $ = id => document.getElementById(id);
const clamp = (v,a,b) => Math.max(a,Math.min(b,v));
const fmt = (v,n=1) => Number.isFinite(v) ? v.toFixed(n) : '—';
const state = {scenario:'elevation',model:'spherical',range:55,nugget:0,sill:400,design:'random',count:8,seed:7,probe:{x:.57,y:.48},samples:[],mode:null,pair:[],drag:-1,grid:null,inverse:null};
const maps = ['krigMap','sdMap','idwMap','truthMap'];
const mapRes = 54;
const rampValue = [[39,79,135],[74,154,188],[184,214,187],[241,205,122],[198,76,55]];
const rampSd = [[234,244,240],[148,197,188],[54,129,140],[41,61,117]];
function rng(seed){let a=seed>>>0;return()=>{a+=0x6D2B79F5;let t=a;t=Math.imul(t^(t>>>15),t|1);t^=t+Math.imul(t^(t>>>7),t|61);return ((t^(t>>>14))>>>0)/4294967296}}
function truth(x,y){
  let v;
  if(state.scenario==='elevation') v=18+49*(.68*x+.42*y)+18*Math.sin(Math.PI*x)*Math.sin(1.5*Math.PI*y)+8*Math.exp(-((x-.73)**2+(y-.72)**2)/.018);
  else if(state.scenario==='temperature') v=46+14*Math.sin(1.4*Math.PI*x+.2)+17*Math.cos(1.2*Math.PI*y-.7*x)+9*Math.sin(3*Math.PI*x+2*Math.PI*y);
  else v=20+67*Math.exp(-((x-.72)**2+(y-.42)**2)/.022)+13*y+7*Math.sin(2.3*Math.PI*x);
  return clamp(v,0,100);
}
function makeSamples(){
  const random=rng(state.seed*9176+(state.scenario==='hotspot'?111:state.scenario==='temperature'?222:333));
  const positions=[];
  if(state.design==='systematic'){
    const cols=4,rows=state.count/4;
    for(let j=0;j<rows;j++) for(let i=0;i<cols;i++) positions.push({x:(i+.5)/cols,y:(j+.5)/rows});
  } else if(state.design==='random'){
    let attempts=0;
    while(positions.length<state.count && attempts<2000){
      attempts++;
      const p={x:.08+.84*random(),y:.08+.84*random()};
      if(positions.every(q=>Math.hypot(q.x-p.x,q.y-p.y)>.11)) positions.push(p);
    }
  } else {
    const centers=[{x:.23,y:.26},{x:.38,y:.74}];
    for(let i=0;i<state.count;i++){
      const c=centers[i%2];let candidate,attempts=0;
      do{
        const angle=2*Math.PI*random(),radius=.035+.16*Math.sqrt(random());
        candidate={x:clamp(c.x+radius*Math.cos(angle),.04,.96),y:clamp(c.y+radius*Math.sin(angle),.04,.96)};
        attempts++;
      }while(attempts<250 && positions.some(q=>Math.hypot(q.x-candidate.x,q.y-candidate.y)<.058));
      positions.push(candidate);
    }
  }
  state.samples=positions.map((p,i)=>({id:i+1,x:p.x,y:p.y,z:Math.round(truth(p.x,p.y)*10)/10}));
  state.pair=[];
}
function dist(a,b){return 100*Math.hypot(a.x-b.x,a.y-b.y)}
function corr(h){
  const t=h/state.range;
  if(state.model==='spherical') return t>=1?0:1-1.5*t+.5*t*t*t;
  if(state.model==='exponential') return Math.exp(-3*t);
  return Math.exp(-3*t*t);
}
function modelGamma(h){return h===0?0:state.nugget+(state.sill-state.nugget)*(1-corr(h))}
function invert(matrix){
  const n=matrix.length,a=matrix.map((row,i)=>row.concat(Array.from({length:n},(_,j)=>i===j?1:0)));
  for(let k=0;k<n;k++){
    let pivot=k;
    for(let i=k+1;i<n;i++) if(Math.abs(a[i][k])>Math.abs(a[pivot][k])) pivot=i;
    if(Math.abs(a[pivot][k])<1e-12) return null;
    [a[k],a[pivot]]=[a[pivot],a[k]];
    const d=a[k][k];
    for(let j=0;j<2*n;j++) a[k][j]/=d;
    for(let i=0;i<n;i++) if(i!==k){
      const factor=a[i][k];
      if(factor!==0) for(let j=0;j<2*n;j++) a[i][j]-=factor*a[k][j];
    }
  }
  return a.map(row=>row.slice(n));
}
function buildInverse(){
  const n=state.samples.length,partial=state.sill-state.nugget;
  for(const ridge of [1e-8,1e-5,.001]){
    const a=Array.from({length:n+1},()=>Array(n+1).fill(0));
    for(let i=0;i<n;i++){
      for(let j=0;j<n;j++) a[i][j]=partial*corr(dist(state.samples[i],state.samples[j]))+(i===j?state.nugget+ridge:0);
      a[i][n]=1;a[n][i]=1;
    }
    const inv=invert(a);
    if(inv){state.inverse=inv;return}
  }
  state.inverse=null;
}
function predict(p,includeWeights=false){
  const s=state.samples,n=s.length,partial=state.sill-state.nugget;
  const ds=s.map(q=>dist(p,q));
  let idw=Array(n).fill(0);
  const exact=ds.findIndex(d=>d<1e-8);
  if(exact>=0) idw[exact]=1;
  else {let total=0;for(let i=0;i<n;i++){idw[i]=1/Math.max(ds[i],.01)**2;total+=idw[i]}idw=idw.map(w=>w/total)}
  const c=s.map(q=>partial*corr(dist(p,q)));
  const b=c.concat(1);
  const solution=state.inverse?state.inverse.map(row=>row.reduce((sum,v,j)=>sum+v*b[j],0)):Array(n).fill(1/n).concat(0);
  const w=solution.slice(0,n);
  let kv=0,iv=0,cross=0;
  for(let i=0;i<n;i++){kv+=w[i]*s[i].z;iv+=idw[i]*s[i].z;cross+=w[i]*c[i]}
  const variance=Math.max(0,partial-cross-solution[n]);
  const out={krig:kv,idw:iv,sd:Math.sqrt(variance)};
  if(includeWeights){out.krigWeights=w;out.idwWeights=idw;out.distances=ds}
  return out;
}
function computeGrid(){
  buildInverse();
  const total=mapRes*mapRes,grid={krig:new Float32Array(total),idw:new Float32Array(total),sd:new Float32Array(total),truth:new Float32Array(total)};
  for(let j=0;j<mapRes;j++) for(let i=0;i<mapRes;i++){
    const p={x:(i+.5)/mapRes,y:(j+.5)/mapRes},k=j*mapRes+i,q=predict(p);
    grid.krig[k]=q.krig;grid.idw[k]=q.idw;grid.sd[k]=q.sd;grid.truth[k]=truth(p.x,p.y);
  }
  state.grid=grid;
}
function color(v,scale,stops){
  const t=clamp(v/scale,0,1)*(stops.length-1),i=Math.min(stops.length-2,Math.floor(t)),f=t-i;
  return [0,1,2].map(k=>Math.round(stops[i][k]*(1-f)+stops[i+1][k]*f));
}
function drawBase(canvas,values,scale,stops){
  const rect=canvas.getBoundingClientRect();if(rect.width<2)return null;
  const dpr=Math.min(window.devicePixelRatio||1,2),w=Math.max(1,Math.round(rect.width*dpr)),h=Math.max(1,Math.round(rect.height*dpr));
  if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h}
  const ctx=canvas.getContext('2d');
  const tiny=document.createElement('canvas');tiny.width=mapRes;tiny.height=mapRes;
  const tc=tiny.getContext('2d'),pixels=tc.createImageData(mapRes,mapRes);
  for(let i=0;i<values.length;i++){const c=color(values[i],scale,stops);pixels.data[i*4]=c[0];pixels.data[i*4+1]=c[1];pixels.data[i*4+2]=c[2];pixels.data[i*4+3]=255}
  tc.putImageData(pixels,0,0);ctx.imageSmoothingEnabled=true;ctx.clearRect(0,0,w,h);ctx.drawImage(tiny,0,0,w,h);
  return {ctx,w,h};
}
function drawMap(id,values,scale,stops){
  const canvas=$(id),base=drawBase(canvas,values,scale,stops);if(!base)return;
  const {ctx,w,h}=base;
  if($('showLinks').checked && id!=='truthMap'){
    ctx.strokeStyle='rgba(15,48,61,.25)';ctx.lineWidth=Math.max(1,w/380);
    for(const s of state.samples){ctx.beginPath();ctx.moveTo(s.x*w,s.y*h);ctx.lineTo(state.probe.x*w,state.probe.y*h);ctx.stroke()}
  }
  if(state.pair.length===2){
    const a=state.samples[state.pair[0]],b=state.samples[state.pair[1]];
    if(a&&b){ctx.strokeStyle='#a6359a';ctx.lineWidth=Math.max(2,w/150);ctx.beginPath();ctx.moveTo(a.x*w,a.y*h);ctx.lineTo(b.x*w,b.y*h);ctx.stroke()}
  }
  if($('showSamples').checked){
    for(let i=0;i<state.samples.length;i++){
      const s=state.samples[i],x=s.x*w,y=s.y*h,r=clamp(w/35,7,13);
      ctx.beginPath();ctx.arc(x,y,r,0,2*Math.PI);ctx.fillStyle=state.pair.includes(i)?'#a6359a':'#153f4a';ctx.fill();ctx.lineWidth=2;ctx.strokeStyle='white';ctx.stroke();
      ctx.font='700 '+Math.round(clamp(w/46,8,12))+'px system-ui';ctx.textAlign='center';ctx.textBaseline='middle';ctx.fillStyle='white';ctx.fillText(String(s.id),x,y+.3);
    }
  }
  const x=state.probe.x*w,y=state.probe.y*h;
  ctx.strokeStyle='white';ctx.lineWidth=5;ctx.beginPath();ctx.arc(x,y,Math.max(7,w/60),0,2*Math.PI);ctx.stroke();
  ctx.strokeStyle='#152b42';ctx.lineWidth=2;ctx.beginPath();ctx.arc(x,y,Math.max(7,w/60),0,2*Math.PI);ctx.stroke();
  ctx.beginPath();ctx.moveTo(x-13,y);ctx.lineTo(x+13,y);ctx.moveTo(x,y-13);ctx.lineTo(x,y+13);ctx.stroke();
}
function drawMaps(){
  if(!state.grid)return;
  drawMap('krigMap',state.grid.krig,100,rampValue);
  drawMap('sdMap',state.grid.sd,30,rampSd);
  drawMap('idwMap',state.grid.idw,100,rampValue);
  if($('showTruth').checked) drawMap('truthMap',state.grid.truth,100,rampValue);
}
function renderReadouts(){
  const p=predict(state.probe,true);
  $('probeCoord').textContent=fmt(state.probe.x*100,0)+', '+fmt(state.probe.y*100,0)+' m';
  $('krigValue').textContent=fmt(p.krig);
  $('sdValue').textContent=fmt(p.sd);
  $('idwValue').textContent=fmt(p.idw);
  $('krigNote').textContent=state.nugget===0?'Exact at observations':'Approximate with noisy observations';
  $('sampleBadge').textContent=state.samples.length+' observations';
  renderWeights(p);
}
function barHtml(weight,signed,maxAbs){
  const percent=Math.min(49,Math.abs(weight)/maxAbs*46);
  const left=signed?(weight<0?50-percent:50):0;
  const width=signed?percent:clamp(weight*100,0,100);
  const cls=signed?(weight<0?'krig-neg':'krig-pos'):'';
  return '<div class="bar-cell"><div class="bar-track '+(signed?'krig':'')+'"><i class="bar-fill '+cls+'" style="left:'+left+'%;width:'+width+'%"></i></div><output>'+fmt(weight,3)+'</output></div>';
}
function renderWeights(p){
  const maxAbs=Math.max(.5,...p.krigWeights.map(Math.abs));
  $('weightBody').innerHTML=state.samples.map((s,i)=>'<tr><td>#'+String(s.id).padStart(2,'0')+'</td><td>'+fmt(s.z)+'</td><td>'+fmt(p.distances[i],0)+' m</td><td>'+barHtml(p.idwWeights[i],false,1)+'</td><td>'+barHtml(p.krigWeights[i],true,maxAbs)+'</td></tr>').join('');
}
function allPairs(){
  const out=[];
  for(let i=0;i<state.samples.length;i++)for(let j=i+1;j<state.samples.length;j++){
    const a=state.samples[i],b=state.samples[j];
    out.push({i,j,h:dist(a,b),gamma:.5*(a.z-b.z)**2});
  }
  return out;
}
function renderPair(pairs){
  if(state.pair.length<2){$('pairValues').textContent=state.pair.length===1?'Sample #'+state.samples[state.pair[0]].id+' selected. Click a second sample.':'Select a pair on a map.';$('pairMath').textContent='γ(h) = ½(zᵢ − zⱼ)²';$('pairDistance').textContent='Distance h and semivariance will appear here.';return null}
  const a=state.samples[state.pair[0]],b=state.samples[state.pair[1]];
  if(!a||!b)return null;
  const chosen=pairs.find(p=>(p.i===state.pair[0]&&p.j===state.pair[1])||(p.j===state.pair[0]&&p.i===state.pair[1]));
  $('pairValues').textContent='Sample #'+a.id+' = '+fmt(a.z)+' · Sample #'+b.id+' = '+fmt(b.z);
  $('pairMath').textContent='γ(h) = ½('+fmt(a.z)+' − '+fmt(b.z)+')² = '+fmt(chosen.gamma);
  $('pairDistance').textContent='Separation h = '+fmt(chosen.h,0)+' m';
  return chosen;
}
function drawVariogram(){
  const canvas=$('variogram'),rect=canvas.getBoundingClientRect();if(rect.width<2)return;
  const dpr=Math.min(window.devicePixelRatio||1,2),w=Math.round(rect.width*dpr),h=Math.round(rect.height*dpr);
  if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h}
  const ctx=canvas.getContext('2d');ctx.clearRect(0,0,w,h);
  const pairs=allPairs(),selected=renderPair(pairs),bins=Array.from({length:10},()=>({sum:0,n:0}));
  for(const p of pairs){const index=Math.min(9,Math.floor(p.h/15));bins[index].sum+=p.gamma;bins[index].n++}
  const peak=Math.max(state.sill*1.2,selected?selected.gamma*1.1:0,...bins.map(b=>b.n?b.sum/b.n*1.1:0),100);
  const yMax=Math.ceil(peak/100)*100;
  const m={left:48*dpr,right:15*dpr,top:17*dpr,bottom:38*dpr},pw=w-m.left-m.right,ph=h-m.top-m.bottom;
  const X=v=>m.left+v/150*pw,Y=v=>m.top+(1-v/yMax)*ph;
  ctx.strokeStyle='#dbe6de';ctx.lineWidth=dpr;
  ctx.fillStyle='#62786f';ctx.font=(10*dpr)+'px system-ui';ctx.textAlign='right';ctx.textBaseline='middle';
  for(let k=0;k<=4;k++){const v=k*yMax/4;ctx.beginPath();ctx.moveTo(m.left,Y(v));ctx.lineTo(w-m.right,Y(v));ctx.stroke();ctx.fillText(String(Math.round(v)),m.left-6*dpr,Y(v))}
  ctx.textAlign='center';ctx.textBaseline='top';
  for(const v of [0,50,100,150]){ctx.fillText(String(v),X(v),h-m.bottom+7*dpr)}
  ctx.fillText('Lag distance (m)',m.left+pw/2,h-14*dpr);
  ctx.save();ctx.translate(13*dpr,m.top+ph/2);ctx.rotate(-Math.PI/2);ctx.textAlign='center';ctx.fillText('Semivariance',0,0);ctx.restore();
  ctx.strokeStyle='#789088';ctx.lineWidth=1.3*dpr;ctx.beginPath();ctx.moveTo(m.left,m.top);ctx.lineTo(m.left,h-m.bottom);ctx.lineTo(w-m.right,h-m.bottom);ctx.stroke();
  if($('showPairs').checked){
    ctx.fillStyle='rgba(53,126,156,.22)';
    for(const p of pairs){ctx.beginPath();ctx.arc(X(p.h),Y(p.gamma),2.7*dpr,0,2*Math.PI);ctx.fill()}
  }
  ctx.strokeStyle='#c56645';ctx.lineWidth=2.5*dpr;ctx.beginPath();
  for(let k=0;k<=150;k++){const x=X(k),y=Y(modelGamma(k));if(k===0)ctx.moveTo(x,y);else ctx.lineTo(x,y)}
  ctx.stroke();
  ctx.fillStyle='#357e9c';
  for(let k=0;k<bins.length;k++)if(bins[k].n){
    const x=X((k+.5)*15),y=Y(bins[k].sum/bins[k].n);
    ctx.beginPath();ctx.arc(x,y,4.5*dpr,0,2*Math.PI);ctx.fill();ctx.strokeStyle='white';ctx.lineWidth=1.5*dpr;ctx.stroke();
  }
  if(selected){
    const x=X(selected.h),y=Y(selected.gamma);
    ctx.beginPath();ctx.arc(x,y,7*dpr,0,2*Math.PI);ctx.fillStyle='#a6359a';ctx.fill();ctx.strokeStyle='white';ctx.lineWidth=2*dpr;ctx.stroke();
  }
  const rangeX=X(state.range);
  ctx.setLineDash([4*dpr,4*dpr]);ctx.strokeStyle='#7d9887';ctx.lineWidth=dpr;ctx.beginPath();ctx.moveTo(rangeX,m.top);ctx.lineTo(rangeX,h-m.bottom);ctx.stroke();ctx.setLineDash([]);
  ctx.font=(10*dpr)+'px system-ui';ctx.fillStyle='#477666';ctx.textAlign='left';ctx.textBaseline='top';ctx.fillText('Range '+state.range+' m',Math.min(rangeX+3*dpr,w-95*dpr),m.top+3*dpr);
}
function renderValidation(){
  if(!$('showTruth').checked)return;
  let krigSq=0,idwSq=0,n=0;
  for(let j=0;j<13;j++)for(let i=0;i<13;i++){
    const x=(i+.5)/13,y=(j+.5)/13,p=predict({x,y}),z=truth(x,y);
    krigSq+=(p.krig-z)**2;idwSq+=(p.idw-z)**2;n++;
  }
  $('krigRmse').textContent=fmt(Math.sqrt(krigSq/n),2);
  $('idwRmse').textContent=fmt(Math.sqrt(idwSq/n),2);
}
function render(full=true){
  if(full)computeGrid();
  drawMaps();renderReadouts();drawVariogram();renderValidation();
}
function setMode(mode){
  if(mode==='pair' && !$('showSamples').checked){$('showSamples').checked=true}
  state.mode=state.mode===mode?null:mode;
  if(state.mode==='pair')state.pair=[];
  $('selectPair').setAttribute('aria-pressed',String(state.mode==='pair'));
  $('addSample').setAttribute('aria-pressed',String(state.mode==='add'));
  $('mapStatus').textContent=state.mode==='pair'?'Pair mode: click two numbered observations on any map.':state.mode==='add'?'Add mode: click an empty location on any map.':'Click a map cell to choose a prediction location, or drag a numbered observation.';
  drawMaps();drawVariogram();
}
function pointFromEvent(e){
  const r=e.currentTarget.getBoundingClientRect();
  return {x:clamp((e.clientX-r.left)/r.width,0,1),y:clamp((e.clientY-r.top)/r.height,0,1),width:r.width};
}
function hitSample(p){
  if(!$('showSamples').checked)return -1;
  let best=-1,delta=1e9;
  state.samples.forEach((s,i)=>{const d=Math.hypot((s.x-p.x)*p.width,(s.y-p.y)*p.width);if(d<delta){delta=d;best=i}});
  return delta<=15?best:-1;
}
function onMapDown(e){
  const p=pointFromEvent(e),hit=hitSample(p);
  if(state.mode==='add'){
    if(state.samples.length>=24){$('mapStatus').textContent='Maximum 24 observations. Reset or change the sampling design to start again.';setMode(null);return}
    state.samples.push({id:state.samples.length+1,x:p.x,y:p.y,z:Math.round(truth(p.x,p.y)*10)/10});
    state.pair=[];setMode(null);render();$('mapStatus').textContent='Added sample #'+state.samples.length+' from the synthetic reference field.';return;
  }
  if(state.mode==='pair'){
    if(hit<0){$('mapStatus').textContent='Click a numbered observation to select a pair.';return}
    if(state.pair.length>=2)state.pair=[];
    if(!state.pair.includes(hit))state.pair.push(hit);
    if(state.pair.length===2){setMode(null);$('mapStatus').textContent='Pair selected. Its semivariance is highlighted below.'}
    else $('mapStatus').textContent='First sample selected. Click a second numbered observation.';
    drawMaps();drawVariogram();return;
  }
  if(hit>=0){
    state.probe={x:state.samples[hit].x,y:state.samples[hit].y};
    render(false);
    state.drag=hit;e.currentTarget.setPointerCapture(e.pointerId);
    $('mapStatus').textContent='Sample #'+state.samples[hit].id+' selected. Drag to move it; its measured value stays the same.';
  } else {
    state.probe={x:p.x,y:p.y};render(false);
    $('mapStatus').textContent='Prediction location: '+fmt(p.x*100,0)+', '+fmt(p.y*100,0)+' m.';
  }
}
function onMapMove(e){
  const p=pointFromEvent(e);
  if(state.drag>=0){state.samples[state.drag].x=p.x;state.samples[state.drag].y=p.y;render();return}
  const hit=hitSample(p);
  if(hit>=0 && !state.mode){const s=state.samples[hit];$('mapStatus').textContent='Sample #'+s.id+': value '+fmt(s.z)+' at '+fmt(s.x*100,0)+', '+fmt(s.y*100,0)+' m. Drag to move it.'}
}
function onMapUp(){state.drag=-1}
function reset(){
  state.scenario='elevation';state.model='spherical';state.range=55;state.nugget=0;state.sill=400;state.design='random';state.count=8;state.seed=7;state.probe={x:.57,y:.48};state.mode=null;state.pair=[];state.drag=-1;
  for(const id of ['scenario','model','design'])$(id).value=state[id];
  for(const id of ['range','nugget','sill','count'])$(id).value=state[id];
  $('showTruth').checked=false;$('showPairs').checked=false;$('showLinks').checked=true;$('showSamples').checked=true;
  $('truthCard').hidden=true;$('validationPanel').hidden=true;$('selectPair').setAttribute('aria-pressed','false');$('addSample').setAttribute('aria-pressed','false');
  $('rangeValue').textContent='55 m';$('nuggetValue').textContent='0';$('sillValue').textContent='400';$('countValue').textContent='8';
  $('mapStatus').textContent='Click a map cell to choose a prediction location, or drag a numbered observation.';
  makeSamples();render();
}
for(const id of maps){const canvas=$(id);canvas.addEventListener('pointerdown',onMapDown);canvas.addEventListener('pointermove',onMapMove);canvas.addEventListener('pointerup',onMapUp);canvas.addEventListener('pointercancel',onMapUp)}
$('selectPair').addEventListener('click',()=>setMode('pair'));
$('addSample').addEventListener('click',()=>setMode('add'));
$('randomize').addEventListener('click',()=>{state.seed++;makeSamples();render();$('mapStatus').textContent='New '+state.design+' sample layout. The synthetic reference field did not change.'});
$('reset').addEventListener('click',reset);
for(const id of ['range','nugget','sill'])$(id).addEventListener('input',e=>{state[id]=Number(e.target.value);$(id+'Value').textContent=String(state[id])+(id==='range'?' m':'');render()});
$('model').addEventListener('change',e=>{state.model=e.target.value;render()});
for(const id of ['scenario','design','count'])$(id).addEventListener(id==='count'?'input':'change',e=>{state[id]=id==='count'?Number(e.target.value):e.target.value;if(id==='count')$('countValue').textContent=String(state.count);makeSamples();render()});
for(const id of ['showLinks','showSamples'])$(id).addEventListener('change',drawMaps);
$('showPairs').addEventListener('change',drawVariogram);
$('showTruth').addEventListener('change',e=>{$('truthCard').hidden=!e.target.checked;$('validationPanel').hidden=!e.target.checked;requestAnimationFrame(()=>{drawMaps();renderValidation()})});
window.addEventListener('resize',()=>{drawMaps();drawVariogram()});
reset();
})();
