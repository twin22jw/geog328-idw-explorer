(() => {
'use strict';
const $ = id => document.getElementById(id);
const clamp = (v,a,b) => Math.max(a,Math.min(b,v));
const fmt = (v,n=1) => Number.isFinite(v) ? v.toFixed(n) : '—';
const rangeLabel = fit => fit.range===20?'≤20 m (fit limit)':fit.range===140?'≥140 m (fit limit)':fit.range+' m';
const state = {range:55,nugget:10,sill:180,design:'random',count:12,seed:7,probe:{x:.57,y:.48},samples:[],mode:null,pair:[],selectedPairs:[],fit:null,lab2:{samples:[],pair:[],selectedPairs:[],reveal:{pairs:false,bins:false,model:false},binIndex:0,fit:null},drag:-1,grid:null,inverse:null,reveal:{pairs:false,bins:false,model:false},binIndex:0};
const maps = ['krigMap','sdMap','idwMap','pairMap','pairMap2'];
const lab = which => which===2?state.lab2:state;
const el = (id,which=1) => $(id+(which===2?'2':''));
const mapRes = 54;
const rampValue = [[39,79,135],[74,154,188],[184,214,187],[241,205,122],[198,76,55]];
const rampSd = [[234,244,240],[148,197,188],[54,129,140],[41,61,117]];
function rng(seed){let a=seed>>>0;return()=>{a+=0x6D2B79F5;let t=a;t=Math.imul(t^(t>>>15),t|1);t^=t+Math.imul(t^(t>>>7),t|61);return ((t^(t>>>14))>>>0)/4294967296}}
function shortField(x,y){
  return clamp(48+19*Math.sin(5*Math.PI*x)*Math.cos(4*Math.PI*y)+11*Math.cos(6*Math.PI*x+2*Math.PI*y),0,100);
}
function longField(x,y){
  return clamp(48+22*Math.sin(1.2*Math.PI*x)*Math.cos(.8*Math.PI*y)+7*Math.cos(1.1*Math.PI*x),0,100);
}
const exampleA={
  points:[[.1527,.6888],[.7902,.5966],[.6926,.4428],[.6572,.3228],[.6272,.5815],[.3767,.2221],[.6101,.9246],[.2767,.0739],[.3003,.6913],[.3747,.4425],[.5934,.6825],[.4761,.8856]],
  values:[63.8,37.2,47.6,45.9,31.4,30.2,60.3,43.7,55.9,40.1,29.7,62.3]
};
const exampleB={
  points:[[.0963,.4765],[.4733,.6033],[.3206,.3146],[.92,.5733],[.924,.853],[.8218,.2602],[.2148,.3703],[.6508,.4664],[.5381,.496],[.4198,.3501],[.7291,.6587],[.7895,.4249]],
  values:[52.9,51.9,48.3,54.8,62.7,49.3,45.4,45.0,45.4,41.9,44.9,38.9]
};
function makeExamples(){
  for(const [L,example] of [[state,exampleA],[state.lab2,exampleB]]){
    L.samples=example.points.map(([x,y],i)=>({id:i+1,x,y,z:example.values[i]}));
    L.pair=[];L.selectedPairs=[];L.reveal={pairs:false,bins:false,model:false};L.binIndex=0;
    L.fit=fitVariogram(L.samples);
  }
  for(const which of [1,2]){el('lagBin',which).value='0';updateRevealControls(which)}
  $('parameterValues').textContent='Case A: range '+rangeLabel(state.fit)+', sill '+fmt(state.fit.sill,1)+', nugget 10. Case B: range '+rangeLabel(state.lab2.fit)+', sill '+fmt(state.lab2.fit.sill,1)+', nugget 10.';
}
function makeSamples(){
  const random=rng(state.seed*9176+333);
  const positions=[];
  if(state.design==='systematic'){
    const cols=4,rows=state.count/4;
    for(let j=0;j<rows;j++) for(let i=0;i<cols;i++) positions.push({x:(i+.5)/cols,y:(j+.5)/rows});
  } else if(state.design==='random'){
    let attempts=0;
    while(positions.length<state.count && attempts<4000){
      attempts++;const p={x:.07+.86*random(),y:.07+.86*random()};
      if(positions.every(q=>Math.hypot(q.x-p.x,q.y-p.y)>.095)) positions.push(p);
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
  state.mapSamples=positions.map((p,i)=>({id:i+1,x:p.x,y:p.y,z:Math.round(longField(p.x,p.y)*10)/10}));
}
function dist(a,b){return 100*Math.hypot(a.x-b.x,a.y-b.y)}
function corr(h){
  const t=h/state.range;
  return t>=1?0:1-1.5*t+.5*t*t*t;
}
function modelGamma(h,fit){const t=h/fit.range;const c=t>=1?0:1-1.5*t+.5*t*t*t;return h===0?0:fit.nugget+(fit.sill-fit.nugget)*(1-c)}
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
  const n=state.mapSamples.length,partial=state.sill-state.nugget;
  for(const ridge of [1e-8,1e-5,.001]){
    const a=Array.from({length:n+1},()=>Array(n+1).fill(0));
    for(let i=0;i<n;i++){
      for(let j=0;j<n;j++) a[i][j]=partial*corr(dist(state.mapSamples[i],state.mapSamples[j]))+(i===j?state.nugget+ridge:0);
      a[i][n]=1;a[n][i]=1;
    }
    const inv=invert(a);
    if(inv){state.inverse=inv;return}
  }
  state.inverse=null;
}
function predict(p,includeWeights=false){
  const s=state.mapSamples,n=s.length,partial=state.sill-state.nugget;
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
  const variance=Math.max(0,state.sill-cross-solution[n]);
  const out={krig:kv,idw:iv,sd:Math.sqrt(variance)};
  if(includeWeights){out.krigWeights=w;out.idwWeights=idw;out.distances=ds}
  return out;
}
function computeGrid(){
  buildInverse();
  const total=mapRes*mapRes,grid={krig:new Float32Array(total),idw:new Float32Array(total),sd:new Float32Array(total)};
  for(let j=0;j<mapRes;j++) for(let i=0;i<mapRes;i++){
    const p={x:(i+.5)/mapRes,y:(j+.5)/mapRes},k=j*mapRes+i,q=predict(p);
    grid.krig[k]=q.krig;grid.idw[k]=q.idw;grid.sd[k]=q.sd;
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
    for(const s of state.mapSamples){ctx.beginPath();ctx.moveTo(s.x*w,s.y*h);ctx.lineTo(state.probe.x*w,state.probe.y*h);ctx.stroke()}
  }
  if($('showSamples').checked){
    for(let i=0;i<state.mapSamples.length;i++){
      const s=state.mapSamples[i],x=s.x*w,y=s.y*h,r=clamp(w/35,7,13);
      ctx.beginPath();ctx.arc(x,y,r,0,2*Math.PI);ctx.fillStyle='#153f4a';ctx.fill();ctx.lineWidth=2;ctx.strokeStyle='white';ctx.stroke();
      ctx.font='700 '+Math.round(clamp(w/46,8,12))+'px system-ui';ctx.textAlign='center';ctx.textBaseline='middle';ctx.fillStyle='white';ctx.fillText(String(s.id),x,y+.3);
    }
  }
  const x=state.probe.x*w,y=state.probe.y*h;
  ctx.strokeStyle='white';ctx.lineWidth=5;ctx.beginPath();ctx.arc(x,y,Math.max(7,w/60),0,2*Math.PI);ctx.stroke();
  ctx.strokeStyle='#152b42';ctx.lineWidth=2;ctx.beginPath();ctx.arc(x,y,Math.max(7,w/60),0,2*Math.PI);ctx.stroke();
  ctx.beginPath();ctx.moveTo(x-13,y);ctx.lineTo(x+13,y);ctx.moveTo(x,y-13);ctx.lineTo(x,y+13);ctx.stroke();
}
function drawPairMap(which=1){
  const L=lab(which),canvas=el('pairMap',which),rect=canvas.getBoundingClientRect();if(rect.width<2)return;
  const dpr=Math.min(window.devicePixelRatio||1,2),w=Math.round(rect.width*dpr),h=Math.round(rect.height*dpr);
  if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h}
  const ctx=canvas.getContext('2d');ctx.fillStyle='#f5f9f3';ctx.fillRect(0,0,w,h);
  ctx.strokeStyle='#dae9dd';ctx.lineWidth=dpr;
  for(let k=1;k<5;k++){const v=k/5;ctx.beginPath();ctx.moveTo(v*w,0);ctx.lineTo(v*w,h);ctx.moveTo(0,v*h);ctx.lineTo(w,v*h);ctx.stroke()}
  if(L.reveal.bins){
    ctx.strokeStyle='rgba(28,104,145,.30)';ctx.lineWidth=1.7*dpr;
    for(const pair of allPairs(L.samples))if(Math.min(9,Math.floor(pair.h/15))===L.binIndex){
      const a=L.samples[pair.i],b=L.samples[pair.j];ctx.beginPath();ctx.moveTo(a.x*w,a.y*h);ctx.lineTo(b.x*w,b.y*h);ctx.stroke();
    }
  }
  if(L.pair.length===2){
    const a=L.samples[L.pair[0]],b=L.samples[L.pair[1]];
    ctx.strokeStyle='#a6359a';ctx.lineWidth=3*dpr;ctx.beginPath();ctx.moveTo(a.x*w,a.y*h);ctx.lineTo(b.x*w,b.y*h);ctx.stroke();
    ctx.fillStyle='#80337c';ctx.font=(12*dpr)+'px system-ui';ctx.textAlign='center';ctx.fillText(fmt(dist(a,b),0)+' m',(a.x+b.x)*w/2,(a.y+b.y)*h/2-8*dpr);
  }
  for(let i=0;i<L.samples.length;i++){
    const s=L.samples[i],x=s.x*w,y=s.y*h,r=10*dpr,chosen=L.pair.includes(i);
    ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.fillStyle=chosen?'#a6359a':'#153f4a';ctx.fill();ctx.lineWidth=2*dpr;ctx.strokeStyle='white';ctx.stroke();
    ctx.fillStyle='white';ctx.textAlign='center';ctx.textBaseline='middle';ctx.font='800 '+(10*dpr)+'px system-ui';ctx.fillText(String(s.id),x,y+.3*dpr);
    const text=fmt(s.z);ctx.font='700 '+(11*dpr)+'px system-ui';
    const tw=ctx.measureText(text).width+10*dpr,tx=clamp(x+9*dpr,2*dpr,w-tw-2*dpr),ty=clamp(y+10*dpr,2*dpr,h-19*dpr);
    ctx.fillStyle='rgba(255,254,250,.92)';ctx.fillRect(tx,ty,tw,17*dpr);
    ctx.fillStyle='#1c5544';ctx.textAlign='left';ctx.textBaseline='middle';ctx.fillText(text,tx+5*dpr,ty+8.5*dpr);
  }
}
function drawMaps(){
  drawPairMap(1);drawPairMap(2);
  if(!state.grid)return;
  drawMap('krigMap',state.grid.krig,100,rampValue);
  drawMap('sdMap',state.grid.sd,30,rampSd);
  drawMap('idwMap',state.grid.idw,100,rampValue);
}
function renderReadouts(){
  const p=predict(state.probe,true);
  $('probeCoord').textContent=fmt(state.probe.x*100,0)+', '+fmt(state.probe.y*100,0)+' m';
  $('krigValue').textContent=fmt(p.krig);
  $('sdValue').textContent=fmt(p.sd);
  $('idwValue').textContent=fmt(p.idw);
  $('krigNote').textContent=state.nugget===0?'Exact at observations':'Approximate with noisy observations';
  $('sampleBadge').textContent=state.mapSamples.length+' observations';
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
  $('weightBody').innerHTML=state.mapSamples.map((s,i)=>'<tr><td>#'+String(s.id).padStart(2,'0')+'</td><td>'+fmt(s.z)+'</td><td>'+fmt(p.distances[i],0)+' m</td><td>'+barHtml(p.idwWeights[i],false,1)+'</td><td>'+barHtml(p.krigWeights[i],true,maxAbs)+'</td></tr>').join('');
}
function allPairs(samples){
  const out=[];
  for(let i=0;i<samples.length;i++)for(let j=i+1;j<samples.length;j++){
    const a=samples[i],b=samples[j];
    out.push({i,j,key:i+'-'+j,h:dist(a,b),gamma:.5*(a.z-b.z)**2});
  }
  return out;
}
function fitVariogram(samples){
  const bins=Array.from({length:10},()=>({sum:0,n:0}));
  for(const p of allPairs(samples)){const b=bins[Math.min(9,Math.floor(p.h/15))];b.sum+=p.gamma;b.n++}
  const observed=bins.map((b,k)=>({h:(k+.5)*15,y:b.n?b.sum/b.n:0,n:b.n})).filter(b=>b.n);
  let best={loss:Infinity,range:55,nugget:10,sill:120};
  for(let range=20;range<=140;range++){
    const parts=observed.map(b=>{const t=Math.min(1,b.h/range);return {b,f:t>=1?1:1.5*t-.5*t*t*t}});
    const denom=parts.reduce((sum,v)=>sum+v.b.n*v.f*v.f,0);
    if(denom<1e-8)continue;
    const sill=10+clamp(parts.reduce((sum,v)=>sum+v.b.n*(v.b.y-10)*v.f,0)/denom,10,890);
    const loss=parts.reduce((sum,v)=>sum+v.b.n*(v.b.y-(10+(sill-10)*v.f))**2,0);
    if(loss<best.loss)best={loss,range,nugget:10,sill};
  }
  return best;
}
function renderPair(pairs,which=1){
  const L=lab(which),first=L.pair.length===1?L.samples[L.pair[0]]:null;
  if(L.pair.length<2){
    el('pairValues',which).textContent=first?'Sample #'+first.id+' = '+fmt(first.z)+'. Choose a second observation.':'Choose two observations on the map.';
    for(const id of ['pairDistance','pairDifference','pairSquared','pairMath'])el(id,which).textContent='—';
    el('pairMapStatus',which).textContent=first?'Sample #'+first.id+' selected. Click a second point.':'Choose your first observation.';
    return null;
  }
  const a=L.samples[L.pair[0]],b=L.samples[L.pair[1]];
  const chosen=pairs.find(p=>(p.i===L.pair[0]&&p.j===L.pair[1])||(p.j===L.pair[0]&&p.i===L.pair[1]));
  const difference=Math.abs(a.z-b.z),squared=difference*difference;
  el('pairValues',which).textContent='Sample #'+a.id+' = '+fmt(a.z)+' · Sample #'+b.id+' = '+fmt(b.z);
  el('pairDistance',which).textContent=fmt(chosen.h,0)+' m';
  el('pairDifference',which).textContent=fmt(difference);
  el('pairSquared',which).textContent=fmt(squared,2);
  el('pairMath',which).textContent='½ × '+fmt(squared,2)+' = '+fmt(chosen.gamma,2);
  el('pairMapStatus',which).textContent='Pair #'+a.id+'–#'+b.id+' is highlighted. Choose another pair to add a point.';
  return chosen;
}
function updateRevealControls(which=1){
  const L=lab(which);
  for(const [key,id] of [['pairs','revealPairs'],['bins','revealBins'],['model','revealModel']])el(id,which).setAttribute('aria-pressed',String(L.reveal[key]));
  el('revealBins',which).disabled=!L.reveal.pairs;
  el('revealModel',which).disabled=!L.reveal.bins;
  el('lagBin',which).disabled=!L.reveal.bins;
  for(const name of ['Range','Nugget','Sill'])el('show'+name,which).disabled=!L.reveal.model;
  el('modelExplorer',which).hidden=!L.reveal.model;
  el('fitReadout',which).textContent='Range '+rangeLabel(L.fit)+' · sill '+fmt(L.fit.sill,1)+' · nugget '+fmt(L.fit.nugget,0);
  el('variogramLegend',which).innerHTML='<span><i class="legend-dot selected"></i>Pairs you selected</span>'
   +(L.reveal.pairs?'<span><i class="legend-dot all"></i>All observation pairs</span>':'')
   +(L.reveal.bins?'<span><i class="legend-dot bin"></i>Experimental lag-bin averages</span>':'')
   +(L.reveal.model?'<span><i class="legend-line fit"></i>Fitted spherical model</span>':'');
}
function renderSelectedPairs(pairs,which=1){
  const L=lab(which),selected=L.selectedPairs.map(key=>pairs.find(p=>p.key===key)).filter(Boolean);
  el('selectedPairCount',which).textContent=selected.length+(selected.length===1?' pair plotted':' pairs plotted');
  el('clearPairs',which).disabled=selected.length===0;
  el('pairHistory',which).textContent=selected.length?'Your points: '+selected.slice(-10).map(p=>'#'+L.samples[p.i].id+'–#'+L.samples[p.j].id).join(', ')+(selected.length>10?' · showing last 10':'')+'.':'Choose two observations to add your first point.';
  return selected;
}
function renderBin(bins,which=1){
  const L=lab(which);
  if(!L.reveal.bins){
    el('binCount',which).textContent='Reveal the lag bins to inspect a distance class.';
    el('binMean',which).textContent='—';el('binFormula',which).textContent='';el('binPairs',which).textContent='';el('binMapHint',which).textContent='';return;
  }
  const k=L.binIndex,b=bins[k],from=k*15,to=(k+1)*15;
  el('binCount',which).textContent=from+'–'+to+' m · '+b.n+(b.n===1?' pair':' pairs');
  el('binMean',which).textContent=b.n?'Experimental mean semivariance = '+fmt(b.sum/b.n,2):'No pairs in this bin';
  el('binFormula',which).textContent=b.n?'Σ pair semivariances = '+fmt(b.sum,2)+'; '+fmt(b.sum,2)+' ÷ '+b.n+' = '+fmt(b.sum/b.n,2):'';
  el('binPairs',which).innerHTML=b.values.slice(0,8).map(p=>'<div><span>#'+L.samples[p.i].id+'–#'+L.samples[p.j].id+' · '+fmt(p.h,0)+' m</span><strong>'+fmt(p.gamma,2)+'</strong></div>').join('')+(b.n>8?'<small>+'+(b.n-8)+' more pairs in this bin</small>':'');
  el('binMapHint',which).textContent=b.n?'These '+b.n+' pair lines are highlighted in blue on the observation map. Their points are dark blue in the graph.':'Try another distance class to see its contributing pairs.';
}
function drawVariogram(which=1){
  const L=lab(which);
  const canvas=el('variogram',which),rect=canvas.getBoundingClientRect();if(rect.width<2)return;
  const dpr=Math.min(window.devicePixelRatio||1,2),w=Math.round(rect.width*dpr),h=Math.round(rect.height*dpr);
  if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h}
  const ctx=canvas.getContext('2d');ctx.clearRect(0,0,w,h);
  const pairs=allPairs(L.samples),selected=renderPair(pairs,which),chosenPairs=renderSelectedPairs(pairs,which),bins=Array.from({length:10},()=>({sum:0,n:0,values:[]}));
  for(const p of pairs){const index=Math.min(9,Math.floor(p.h/15));bins[index].sum+=p.gamma;bins[index].n++;bins[index].values.push(p)}
  renderBin(bins,which);
  const peak=Math.max(L.fit.sill*1.5,selected?selected.gamma*1.1:0,...bins.map(b=>b.n?b.sum/b.n*1.15:0),...chosenPairs.map(p=>p.gamma*1.1),60);
  const yMax=Math.ceil(peak/25)*25;
  const m={left:50*dpr,right:16*dpr,top:19*dpr,bottom:41*dpr},pw=w-m.left-m.right,ph=h-m.top-m.bottom;
  const X=v=>m.left+v/150*pw,Y=v=>m.top+(1-clamp(v,0,yMax)/yMax)*ph;
  const dot=(p,r,fill)=>{const clipped=p.gamma>yMax;ctx.fillStyle=fill;ctx.beginPath();ctx.arc(X(p.h),Y(p.gamma)+(clipped?4*dpr:0),r,0,2*Math.PI);ctx.fill();if(clipped){ctx.font=(11*dpr)+'px system-ui';ctx.fillText('↑',X(p.h),m.top-1*dpr)}};
  ctx.strokeStyle='#dbe6de';ctx.lineWidth=dpr;
  ctx.fillStyle='#62786f';ctx.font=(10*dpr)+'px system-ui';ctx.textAlign='right';ctx.textBaseline='middle';
  for(let k=0;k<=4;k++){const v=k*yMax/4;ctx.beginPath();ctx.moveTo(m.left,Y(v));ctx.lineTo(w-m.right,Y(v));ctx.stroke();ctx.fillText(String(Math.round(v)),m.left-6*dpr,Y(v))}
  ctx.textAlign='center';ctx.textBaseline='top';
  for(const v of [0,50,100,150])ctx.fillText(String(v),X(v),h-m.bottom+7*dpr);
  ctx.fillText('Separation distance h (m)',m.left+pw/2,h-14*dpr);
  ctx.save();ctx.translate(13*dpr,m.top+ph/2);ctx.rotate(-Math.PI/2);ctx.textAlign='center';ctx.fillText('Semivariance',0,0);ctx.restore();
  if(L.reveal.bins){
    ctx.fillStyle='rgba(52,126,155,.08)';ctx.fillRect(X(L.binIndex*15),m.top,X(15)-X(0),ph);
  }
  ctx.strokeStyle='#789088';ctx.lineWidth=1.3*dpr;ctx.beginPath();ctx.moveTo(m.left,m.top);ctx.lineTo(m.left,h-m.bottom);ctx.lineTo(w-m.right,h-m.bottom);ctx.stroke();
  if(L.reveal.pairs){
    ctx.fillStyle='rgba(53,126,156,.25)';
    for(const p of pairs)dot(p,2.8*dpr,'rgba(53,126,156,.25)')
  }
  if(L.reveal.bins){
    ctx.fillStyle='rgba(24,96,136,.85)';
    for(const p of bins[L.binIndex].values)dot(p,4*dpr,'rgba(24,96,136,.85)')
    ctx.fillStyle='#357e9c';
    for(let k=0;k<bins.length;k++)if(bins[k].n){
      const x=X((k+.5)*15),y=Y(bins[k].sum/bins[k].n);
      ctx.beginPath();ctx.arc(x,y,k===L.binIndex?6*dpr:4.5*dpr,0,2*Math.PI);ctx.fill();ctx.strokeStyle='white';ctx.lineWidth=1.5*dpr;ctx.stroke();
    }
  }
  if(L.reveal.model){
    ctx.strokeStyle='#c56645';ctx.lineWidth=2.7*dpr;ctx.beginPath();
    for(let k=0;k<=150;k++){const x=X(k),y=Y(modelGamma(k,L.fit));if(k===0)ctx.moveTo(x,y);else ctx.lineTo(x,y)}
    ctx.stroke();
    ctx.font=(11*dpr)+'px system-ui';ctx.textAlign='left';ctx.textBaseline='middle';
    if(el('showSill',which).checked){
      ctx.setLineDash([5*dpr,4*dpr]);ctx.strokeStyle='#27699b';ctx.lineWidth=1.5*dpr;ctx.beginPath();ctx.moveTo(m.left,Y(L.fit.sill));ctx.lineTo(w-m.right,Y(L.fit.sill));ctx.stroke();ctx.setLineDash([]);
      ctx.fillStyle='#195b83';ctx.fillText('Sill '+fmt(L.fit.sill,1),Math.max(m.left+5*dpr,w-m.right-92*dpr),Y(L.fit.sill)-9*dpr);
    }
    if(el('showRange',which).checked){
      const rx=X(L.fit.range);ctx.setLineDash([5*dpr,4*dpr]);ctx.strokeStyle='#24775d';ctx.lineWidth=1.5*dpr;ctx.beginPath();ctx.moveTo(rx,m.top);ctx.lineTo(rx,h-m.bottom);ctx.stroke();ctx.setLineDash([]);
      ctx.fillStyle='#17654e';ctx.fillText('Range '+L.fit.range+' m',Math.min(rx+5*dpr,w-100*dpr),m.top+12*dpr);
      ctx.strokeStyle='#24775d';ctx.lineWidth=2*dpr;const yy=h-m.bottom-12*dpr;ctx.beginPath();ctx.moveTo(X(0),yy);ctx.lineTo(rx,yy);ctx.stroke();ctx.fillText('↔',X(L.fit.range/2)-5*dpr,yy-10*dpr);
    }
    if(el('showNugget',which).checked){
      ctx.strokeStyle='#b24b81';ctx.lineWidth=3*dpr;ctx.beginPath();ctx.moveTo(X(0)+5*dpr,Y(0));ctx.lineTo(X(0)+5*dpr,Y(L.fit.nugget));ctx.stroke();
      ctx.fillStyle='#a33371';ctx.beginPath();ctx.arc(X(0),Y(L.fit.nugget),4*dpr,0,2*Math.PI);ctx.fill();
      ctx.fillText('Nugget '+fmt(L.fit.nugget,0),X(0)+11*dpr,Y(L.fit.nugget)-12*dpr);
    }
  }
  for(const p of chosenPairs){
    ctx.beginPath();ctx.arc(X(p.h),Y(p.gamma),4.5*dpr,0,2*Math.PI);ctx.fillStyle='#a6359a';ctx.fill();ctx.strokeStyle='white';ctx.lineWidth=1.2*dpr;ctx.stroke();
  }
  if(selected){
    const x=X(selected.h),y=Y(selected.gamma);
    ctx.beginPath();ctx.arc(x,y,7*dpr,0,2*Math.PI);ctx.fillStyle='#a6359a';ctx.fill();ctx.strokeStyle='white';ctx.lineWidth=2*dpr;ctx.stroke();
  }else if(!L.reveal.pairs && !chosenPairs.length){
    ctx.fillStyle='#6a8174';ctx.font=(13*dpr)+'px system-ui';ctx.textAlign='center';ctx.textBaseline='middle';
    ctx.fillText('Choose two observations to place the first point',m.left+pw/2,m.top+ph/2);
  }
}
function render(full=true){
  if(full)computeGrid();
  drawMaps();renderReadouts();drawVariogram(1);drawVariogram(2);
}
function pointFromEvent(e){
  const r=e.currentTarget.getBoundingClientRect();
  return {x:clamp((e.clientX-r.left)/r.width,0,1),y:clamp((e.clientY-r.top)/r.height,0,1),width:r.width};
}
function hitSample(p,samples,force=false){
  if(!force && !$('showSamples').checked)return -1;
  let best=-1,delta=1e9;
  samples.forEach((s,i)=>{const d=Math.hypot((s.x-p.x)*p.width,(s.y-p.y)*p.width);if(d<delta){delta=d;best=i}});
  return delta<=15?best:-1;
}
function choosePair(hit,which){
  const L=lab(which);
  if(hit<0){el('pairMapStatus',which).textContent='Click a numbered observation to choose a pair.';return}
  if(L.pair.length>=2)L.pair=[];
  if(!L.pair.includes(hit))L.pair.push(hit);
  if(L.pair.length===2){
    const key=[...L.pair].sort((a,b)=>a-b).join('-');
    if(!L.selectedPairs.includes(key))L.selectedPairs.push(key);
  }
  drawPairMap(which);drawVariogram(which);
}
function onMapDown(e){
  const p=pointFromEvent(e),id=e.currentTarget.id,which=id==='pairMap2'?2:1;
  if(id==='pairMap'||id==='pairMap2'){choosePair(hitSample(p,lab(which).samples,true),which);return}
  const hit=hitSample(p,state.mapSamples);
  if(hit>=0){
    state.probe={x:state.mapSamples[hit].x,y:state.mapSamples[hit].y};
    render(false);
    state.drag=hit;e.currentTarget.setPointerCapture(e.pointerId);
    $('mapStatus').textContent='Sample #'+state.mapSamples[hit].id+' selected. Drag to move it; its value stays fixed in this example.';
  }else{
    state.probe={x:p.x,y:p.y};render(false);
    $('mapStatus').textContent='Prediction location: '+fmt(p.x*100,0)+', '+fmt(p.y*100,0)+' m.';
  }
}
function onMapMove(e){
  const p=pointFromEvent(e),id=e.currentTarget.id,which=id==='pairMap2'?2:1;
  if(state.drag>=0){
    state.mapSamples[state.drag].x=p.x;state.mapSamples[state.drag].y=p.y;render();return;
  }
  const L=(id==='pairMap'||id==='pairMap2')?lab(which):{samples:state.mapSamples},hit=hitSample(p,L.samples,id==='pairMap'||id==='pairMap2');
  if(hit>=0){
    const s=L.samples[hit];
    if(id==='pairMap'||id==='pairMap2')el('pairMapStatus',which).textContent='Sample #'+s.id+' has value '+fmt(s.z)+'. Click to select it.';
    else $('mapStatus').textContent='Sample #'+s.id+': value '+fmt(s.z)+' at '+fmt(s.x*100,0)+', '+fmt(s.y*100,0)+' m.';
  }
}
function onMapUp(){state.drag=-1}
function reset(){
  state.range=55;state.nugget=10;state.sill=180;state.design='random';state.count=12;state.seed=7;state.probe={x:.57,y:.48};state.drag=-1;
  $('design').value=state.design;$('count').value=String(state.count);$('countValue').textContent='12';
  $('showLinks').checked=true;$('showSamples').checked=true;
  $('mapStatus').textContent='Click a map cell to choose a prediction location, or drag a numbered observation.';
  $('range').value='55';$('rangeValue').textContent='55 m';$('nugget').value='10';$('nuggetValue').textContent='10';$('sill').value='180';$('sillValue').textContent='180';
  for(const which of [1,2])for(const name of ['Range','Nugget','Sill'])el('show'+name,which).checked=false;
  makeExamples();makeSamples();render();
}
for(const id of maps){
  const canvas=$(id);canvas.addEventListener('pointerdown',onMapDown);canvas.addEventListener('pointermove',onMapMove);
  canvas.addEventListener('pointerup',onMapUp);canvas.addEventListener('pointercancel',onMapUp);
}
$('randomize').addEventListener('click',()=>{state.seed++;makeSamples();render();$('mapStatus').textContent='New independent sample layout for the prediction maps.'});
$('reset').addEventListener('click',reset);
$('range').addEventListener('input',e=>{state.range=Number(e.target.value);$('rangeValue').textContent=state.range+' m';render()});
$('nugget').addEventListener('input',e=>{state.nugget=Number(e.target.value);state.sill=Math.max(state.sill,state.nugget+10);$('nuggetValue').textContent=String(state.nugget);$('sill').value=String(state.sill);$('sillValue').textContent=String(state.sill);render()});
$('sill').addEventListener('input',e=>{state.sill=Math.max(Number(e.target.value),state.nugget+10);e.target.value=String(state.sill);$('sillValue').textContent=String(state.sill);render()});
$('design').addEventListener('change',e=>{state.design=e.target.value;makeSamples();render()});
$('count').addEventListener('input',e=>{state.count=Number(e.target.value);$('countValue').textContent=String(state.count);makeSamples();render()});
for(const id of ['showLinks','showSamples'])$(id).addEventListener('change',drawMaps);
for(const which of [1,2]){
  el('clearPairs',which).addEventListener('click',()=>{const L=lab(which);L.selectedPairs=[];L.pair=[];drawPairMap(which);drawVariogram(which)});
  el('revealPairs',which).addEventListener('click',()=>{const L=lab(which);L.reveal.pairs=!L.reveal.pairs;if(!L.reveal.pairs){L.reveal.bins=false;L.reveal.model=false}updateRevealControls(which);drawPairMap(which);drawVariogram(which)});
  el('revealBins',which).addEventListener('click',()=>{const L=lab(which);L.reveal.bins=!L.reveal.bins;if(!L.reveal.bins)L.reveal.model=false;updateRevealControls(which);drawPairMap(which);drawVariogram(which)});
  el('revealModel',which).addEventListener('click',()=>{const L=lab(which);L.reveal.model=!L.reveal.model;updateRevealControls(which);drawVariogram(which)});
  el('lagBin',which).addEventListener('change',e=>{lab(which).binIndex=Number(e.target.value);drawPairMap(which);drawVariogram(which)});
  for(const name of ['Range','Nugget','Sill'])el('show'+name,which).addEventListener('change',()=>drawVariogram(which));
}
$('explainParams').addEventListener('click',()=>$('parameterDialog').showModal());
$('closeParameterDialog').addEventListener('click',()=>$('parameterDialog').close());
$('parameterDialog').addEventListener('click',e=>{if(e.target===$('parameterDialog'))e.target.close()});
window.addEventListener('resize',()=>{drawMaps();drawVariogram(1);drawVariogram(2)});
reset();
})();



