const $ = id => document.getElementById(id);
const colors = ['#247466', '#be714b', '#6579a3', '#a18b3b', '#8d6090', '#338c9d', '#9c624c', '#697d4c'];
const escapeHTML = value => String(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const displayPopulation = name => name.replace(/^--\s*/, '').replace('England & Wales Total Population', 'England & Wales').replace('France total population', 'France').replace('NZ total population', 'New Zealand');
const titleCase = word => word.charAt(0).toUpperCase() + word.slice(1);
const formatRate = value => value.toLocaleString('en-US', {maximumSignificantDigits: 5});
const defaults = () => ({method:'DSA', mode:'age', populations:['Canada'], sexes:['female','male'], age:60, ages:[40,60,80], year:2025, log:false});
let state = defaults();
let manifest, data, plotted = [], records = [], plotGeometry, activeRequest = 0;
const cache = new Map();
const currentMethod = () => manifest.methods.find(m => m.id === state.method);
function readLocation() {
  if (!location.hash.startsWith('#view=')) return;
  try {
    const candidate = JSON.parse(decodeURIComponent(location.hash.slice(6)));
    if (candidate && typeof candidate === 'object') state = {...defaults(), ...candidate};
  } catch { state = defaults(); }
}
function normalize() {
  if (!manifest.methods.some(m => m.id === state.method)) state.method = 'DSA';
  if (!['age','population','year'].includes(state.mode)) state.mode = 'age';
  const m = currentMethod();
  const initialPopulations = Array.isArray(state.populations) ? state.populations : [];
  state.populations = [...new Set(initialPopulations.filter(p => m.populations.includes(p)))];
  if (initialPopulations.length && !state.populations.length) state.populations = [m.populations.includes('Canada') ? 'Canada' : m.populations[0]];
  state.sexes = [...new Set((Array.isArray(state.sexes) ? state.sexes : []).filter(s => m.sexes.includes(s)))];
  if (state.mode === 'population') {state.populations = state.populations.slice(0,1); state.sexes = state.sexes.slice(0,1);}
  state.age = m.ages.includes(Number(state.age)) ? Number(state.age) : m.ages.reduce((a,b) => Math.abs(b-Number(state.age || 60)) < Math.abs(a-Number(state.age || 60)) ? b : a);
  state.year = m.years.includes(Number(state.year)) ? Number(state.year) : m.years.at(-1);
  state.ages = [...new Set((Array.isArray(state.ages) ? state.ages : []).map(Number).filter(a => m.ages.includes(a)))].sort((a,b) => a-b);
  state.log = state.log === true;
}
function saveLocation() { history.replaceState(null, '', '#view=' + encodeURIComponent(JSON.stringify(state))); }
function populateControls() {
  const m = currentMethod();
  $('method').value = state.method;
  $('method-hint').textContent = m.description;
  $('coverage-note').textContent = `${m.populations.length} populations · ages ${m.ages[0]}–${m.ages.at(-1)} · ${m.years[0]}–${m.years.at(-1)}. ${state.method === 'CBD' ? 'CBD results are available for older ages only.' : 'The available populations may differ between methods.'}`;
  $('country-hint').textContent = state.mode === 'population' ? 'Choose one population and one sex.' : `${state.populations.length} selected · compare one or more populations.`;
  for (const sex of ['female','male']) {
    $(sex).type = state.mode === 'population' ? 'radio' : 'checkbox';
    $(sex).name = state.mode === 'population' ? 'selected-sex' : sex;
    $(sex).checked = state.sexes.includes(sex);
  }
  $('age').min = m.ages[0]; $('age').max = m.ages.at(-1); $('age').value = state.age;
  $('age-value').value = state.age; $('age-min').textContent = `${m.ages[0]} years`; $('age-max').textContent = `${m.ages.at(-1)} years`;
  $('year').min = m.years[0]; $('year').max = m.years.at(-1); $('year').value = state.year; $('year-value').value = state.year;
  $('ages').value = state.ages.join(', '); $('log-scale').checked = state.log;
  $('ages-hint').textContent = `Use comma-separated ages from ${m.ages[0]} to ${m.ages.at(-1)}.`;
  $('age-field').hidden = state.mode !== 'age'; $('ages-field').hidden = state.mode !== 'population'; $('year-field').hidden = state.mode !== 'year';
  document.querySelectorAll('[data-mode]').forEach(button => {const active = button.dataset.mode === state.mode; button.classList.toggle('active',active); button.setAttribute('aria-pressed', String(active));});
  renderPopulations();
}
function renderPopulations() {
  const query = $('country-search').value.trim().toLowerCase();
  const populations = currentMethod().populations.filter(p => displayPopulation(p).toLowerCase().includes(query)).sort((a,b) => Number(state.populations.includes(b))-Number(state.populations.includes(a)) || displayPopulation(a).localeCompare(displayPopulation(b)));
  $('country-options').replaceChildren();
  for (const p of populations) {
    const label = document.createElement('label'); label.className = 'country-option';
    const input = document.createElement('input'); input.type = state.mode === 'population' ? 'radio' : 'checkbox'; input.name = 'population'; input.checked = state.populations.includes(p); input.value = p;
    input.addEventListener('change', () => {
      state.populations = state.mode === 'population' ? [p] : input.checked ? [...state.populations,p] : state.populations.filter(x=>x!==p);
      $('country-hint').textContent = state.mode === 'population' ? 'Choose one population and one sex.' : `${state.populations.length} selected · compare one or more populations.`;
      render();
    });
    const text = document.createElement('span'); text.textContent = displayPopulation(p);
    label.append(input,text); $('country-options').append(label);
  }
  if (!populations.length) { const empty = document.createElement('p'); empty.className='field-hint'; empty.textContent='No matching populations.'; $('country-options').append(empty); }
}
function setStatus(message, retry = false) {
  $('status').replaceChildren(document.createTextNode(message)); $('status').hidden = false; $('chart').setAttribute('hidden','');
  $('legend').replaceChildren(); $('tooltip').hidden = true;
  for (const id of ['download','save-chart','show-table']) $(id).disabled = true;
  $('table-panel').hidden = true; $('show-table').setAttribute('aria-expanded','false'); $('show-table').textContent = 'View data';
  $('point-count').textContent = '0'; plotted = []; records = [];
  if (retry) {const button = document.createElement('button'); button.textContent='Try again'; button.onclick=()=>manifest ? loadMethod() : init(); $('status').append(button);}
}
async function fetchJSON(path) {const response = await fetch(path); if (!response.ok) throw new Error(`HTTP ${response.status}`); return response.json();}
async function loadMethod() {
  const request = ++activeRequest; const method = state.method;
  setStatus('Loading saved forecasts…');
  try {
    if (!cache.has(method)) cache.set(method, await fetchJSON(`data/${currentMethod().file}`));
    if (request !== activeRequest) return;
    data = cache.get(method); render();
  } catch {if (request === activeRequest) setStatus('The forecast data could not be loaded. Check your connection and try again.', true);}
}
function buildSeries() {
  const list = [];
  for (const population of state.populations) for (const sex of state.sexes) {
    const source = data.series.find(s=>s.population===population && s.sex===sex);
    if (!source) continue;
    if (state.mode === 'year') {
      const yearIndex = data.years.indexOf(state.year);
      list.push({label:`${displayPopulation(population)} · ${titleCase(sex)}`, population, sex, points:data.ages.map((age,i)=>({x:age,age,year:state.year,rate:source.values[i][yearIndex]}))});
    } else {
      const ages = state.mode === 'age' ? [state.age] : state.ages;
      for (const age of ages) {
        const row = source.values[data.ages.indexOf(age)];
        if (!row) continue;
        list.push({label:state.mode === 'population' ? `Age ${age}` : `${displayPopulation(population)} · ${titleCase(sex)}`,population,sex,points:data.years.map((year,i)=>({x:year,age,year,rate:row[i]}))});
      }
    }
  }
  return list.map((s,i)=>({...s,color:colors[i%colors.length],dash:state.mode!=='population' && s.sex==='male' ? '7 4' : i>=colors.length ? '3 3' : ''}));
}
function render() {
  saveLocation(); $('notice').textContent = ''; $('tooltip').hidden = true;
  if (!data || data.method !== state.method) return;
  plotted = buildSeries();
  if (!plotted.length) {setStatus(state.mode==='population' && !state.ages.length ? 'Enter at least one age in the available range to see a forecast.' : 'Choose a population and sex to explore the forecasts.'); return;}
  const subtitle = state.mode === 'population' ? `${displayPopulation(state.populations[0])} · ${titleCase(state.sexes[0])} · ${currentMethod().label}` : `${currentMethod().label} · ${state.populations.length} population${state.populations.length===1?'':'s'} · ${state.mode==='year' ? 'Year '+state.year : 'Age '+state.age}`;
  $('chart-kicker').textContent = state.mode==='year' ? 'AGE PROFILES' : 'FORECAST TRAJECTORIES';
  $('chart-title').textContent = state.mode==='year' ? `Mortality across ages in ${state.year}` : state.mode==='population' ? 'How forecasts vary with age' : `Mortality over time at age ${state.age}`;
  $('chart-subtitle').textContent = subtitle;
  $('selection-summary').textContent = `${plotted.length} ${plotted.length===1?'series':'series'} · ${data.years[0]}–${data.years.at(-1)} forecasts`;
  $('status').hidden = true; $('chart').removeAttribute('hidden');
  for (const id of ['download','save-chart','show-table']) $(id).disabled = false;
  records = plotted.flatMap(s=>s.points.map(p=>({method:state.method,population:s.population,sex:s.sex,age:p.age,year:p.year,mortality_rate:p.rate,rate_times_100:p.rate*100})));
  $('point-count').textContent = records.length.toLocaleString();
  drawChart(); renderTable();
  $('legend').innerHTML = plotted.map(s=>`<span class="legend-item"><span class="legend-swatch" style="border-color:${s.color};border-top-style:${s.dash?'dashed':'solid'}"></span>${escapeHTML(s.label)}</span>`).join('');
}
function niceStep(range) {const rough = range / 5; const power = 10 ** Math.floor(Math.log10(rough)); return [1,2,2.5,5,10].find(n=>n*power>=rough)*power;}
function drawChart() {
  const width = Math.max(390, $('chart-container').clientWidth - 28), height = width < 580 ? 310 : 400;
  const pad = {left:65,right:24,top:38,bottom:50}; const w=width-pad.left-pad.right, h=height-pad.top-pad.bottom;
  const points = plotted.flatMap(s=>s.points), values=points.map(p=>p.rate*100).filter(v=>Number.isFinite(v) && (!state.log || v>0));
  if (!values.length) {setStatus('This selection has no positive rates to display on a log scale.'); return;}
  const min = Math.min(...values), max = Math.max(...values);
  let ymin=0,ymax,ticks=[];
  if (state.log) {
    ymin=Math.floor(Math.log10(min)); ymax=Math.ceil(Math.log10(max)); if(ymin===ymax)ymax++;
    const stride=Math.max(1,Math.ceil((ymax-ymin)/6)); for(let t=ymin;t<=ymax;t+=stride)ticks.push(10**t);
  } else {const step=niceStep(max||1); ymax=Math.ceil(max/step)*step || 1; for(let t=0;t<=ymax+step/100;t+=step)ticks.push(t);}
  const xs=plotted[0].points.map(p=>p.x), xmin=xs[0], xmax=xs.at(-1);
  const x=v=>pad.left+(v-xmin)/(xmax-xmin||1)*w;
  const y=v=>pad.top+h-((state.log?Math.log10(v):v)-ymin)/(ymax-ymin)*h;
  const tickCount = width<580 ? 4 : 7;
  const xticks=[...new Set(Array.from({length:tickCount},(_,i)=>xs[Math.round(i*(xs.length-1)/(tickCount-1))]))];
  const chart=$('chart'); chart.setAttribute('viewBox',`0 0 ${width} ${height}`);
  const axisFont='font-family="Segoe UI,Arial,sans-serif" font-size="11" fill="#78847c"';
  const yLabel=`Mortality rate × 100${state.log?' · log scale':''}`;
  chart.innerHTML=`<title id="svg-title">${escapeHTML($('chart-title').textContent)}</title><desc id="svg-desc">${escapeHTML($('chart-subtitle').textContent)}. ${escapeHTML(yLabel)}. Saved point forecasts. Use left and right arrow keys to inspect values, or select View data for an accessible table.</desc><rect width="${width}" height="${height}" fill="white"/><text x="${pad.left}" y="18" ${axisFont}>${yLabel}</text>`+
    ticks.map(t=>`<line x1="${pad.left}" x2="${width-pad.right}" y1="${y(t)}" y2="${y(t)}" stroke="#e9eee7"/><text x="${pad.left-12}" y="${y(t)+4}" text-anchor="end" ${axisFont}>${formatRate(t)}</text>`).join('')+
    xticks.map(t=>`<text x="${x(t)}" y="${height-pad.bottom+25}" text-anchor="middle" ${axisFont}>${t}</text>`).join('')+
    `<text x="${pad.left+w/2}" y="${height-5}" text-anchor="middle" ${axisFont}>${state.mode==='year'?'Age (years)':'Forecast year'}</text>`+
    plotted.map(s=>{let pen=false; const d=s.points.map(p=>{const v=p.rate*100;if(!Number.isFinite(v)||(state.log&&v<=0)){pen=false;return '';} const command=pen?'L':'M';pen=true;return `${command}${x(p.x).toFixed(2)},${y(v).toFixed(2)}`;}).join(' '); return `<path d="${d}" fill="none" stroke="${s.color}" stroke-width="2.6" stroke-linecap="round" stroke-linejoin="round" ${s.dash?`stroke-dasharray="${s.dash}"`:''}/>`;}).join('')+
    '<g id="crosshair" style="display:none"><line id="crosshair-line" stroke="#8da59c" stroke-dasharray="3 4"/><g id="crosshair-dots"></g></g>';
  plotGeometry={width,height,pad,x,y,xs};
}
let focusIndex = 0;
function inspect(index) {
  if(!plotGeometry||!plotted.length)return;
  const {xs,x,y,pad,width,height}=plotGeometry; focusIndex=Math.max(0,Math.min(index,xs.length-1));
  const line=$('crosshair-line'); if(!line)return;
  $('crosshair').style.display=''; line.setAttribute('x1',x(xs[focusIndex]));line.setAttribute('x2',x(xs[focusIndex]));line.setAttribute('y1',pad.top);line.setAttribute('y2',height-pad.bottom);
  $('crosshair-dots').innerHTML=plotted.filter(s=>!state.log||s.points[focusIndex].rate>0).map(s=>`<circle cx="${x(xs[focusIndex])}" cy="${y(s.points[focusIndex].rate*100)}" r="4" fill="${s.color}" stroke="white" stroke-width="2"/>`).join('');
  const tooltip=$('tooltip'); tooltip.hidden=false;
  tooltip.innerHTML=`<strong>${state.mode==='year'?'Age':'Year'} ${xs[focusIndex]}</strong>`+plotted.slice(0,12).map(s=>`<div class="tooltip-line"><span>${escapeHTML(s.label)}</span><span>${formatRate(s.points[focusIndex].rate*100)}</span></div>`).join('')+(plotted.length>12?'<div>More series in the data table</div>':'');
  const container=$('chart-container'); const offset=x(xs[focusIndex])+22;
  tooltip.style.left=`${Math.max(5,Math.min(offset,container.clientWidth-tooltip.offsetWidth-8))}px`; tooltip.style.top='35px';
}
function hideTooltip(){ $('tooltip').hidden=true; if($('crosshair'))$('crosshair').style.display='none'; }
function renderTable(){if($('table-panel').hidden)return; $('table-body').innerHTML=records.map(r=>`<tr><td>${escapeHTML(displayPopulation(r.population))}</td><td>${titleCase(r.sex)}</td><td>${r.age}</td><td>${r.year}</td><td>${formatRate(r.rate_times_100)}</td></tr>`).join('');}
function downloadBlob(content,type,filename){const url=URL.createObjectURL(new Blob([content],{type}));const a=document.createElement('a');a.href=url;a.download=filename;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function setupEvents(){
  $('method').addEventListener('change',()=>{state.method=$('method').value;normalize();populateControls();saveLocation();loadMethod();});
  $('country-search').addEventListener('input',renderPopulations);
  document.querySelectorAll('[data-mode]').forEach(button=>button.addEventListener('click',()=>{state.mode=button.dataset.mode;normalize();populateControls();render();}));
  for(const sex of ['female','male'])$(sex).addEventListener('change',()=>{state.sexes=['female','male'].filter(s=>$(s).checked);render();});
  $('age').addEventListener('input',()=>{state.age=Number($('age').value);$('age-value').value=state.age;render();});
  $('year').addEventListener('input',()=>{state.year=Number($('year').value);$('year-value').value=state.year;render();});
  $('ages').addEventListener('input',()=>{const raw=$('ages').value.trim(); const requested=raw?raw.split(',').map(s=>s.trim()).filter(Boolean).map(Number):[];state.ages=[...new Set(requested.filter(a=>currentMethod().ages.includes(a)))].sort((a,b)=>a-b);render();if(requested.some(a=>!currentMethod().ages.includes(a)))$('notice').textContent=`Only available integer ages ${data.ages[0]}–${data.ages.at(-1)} are shown.`;});
  $('log-scale').addEventListener('change',()=>{state.log=$('log-scale').checked;render();});
  $('reset').addEventListener('click',()=>{state=defaults();$('country-search').value='';populateControls();saveLocation();loadMethod();});
  $('show-table').addEventListener('click',()=>{const show=$('table-panel').hidden;$('table-panel').hidden=!show;$('show-table').setAttribute('aria-expanded',String(show));$('show-table').textContent=show?'Hide data':'View data';renderTable();});
  $('download').addEventListener('click',()=>{if(!records.length)return;const columns=Object.keys(records[0]);const quote=v=>'"'+String(v).replaceAll('"','""')+'"';const csv=[columns.join(','),...records.map(r=>columns.map(k=>quote(r[k])).join(','))].join('\r\n');downloadBlob(csv,'text/csv;charset=utf-8',`mortselect-${state.method}-${state.mode}.csv`);});
  $('save-chart').addEventListener('click',()=>{if(!plotted.length)return;const copy=$('chart').cloneNode(true);copy.setAttribute('xmlns','http://www.w3.org/2000/svg');copy.removeAttribute('hidden');copy.querySelector('#crosshair')?.remove();const {width,height}=plotGeometry;const rows=Math.ceil(plotted.length/2), totalHeight=height+110+rows*23;copy.setAttribute('viewBox',`0 -54 ${width} ${totalHeight}`);copy.setAttribute('width',width);copy.setAttribute('height',totalHeight);const original=copy.innerHTML;copy.innerHTML=`<rect x="0" y="-54" width="${width}" height="${totalHeight}" fill="white"/><text x="20" y="-29" fill="#183c38" font-family="Arial" font-size="18">${escapeHTML($('chart-title').textContent)}</text><text x="20" y="-8" fill="#697774" font-family="Arial" font-size="11">${escapeHTML($('chart-subtitle').textContent)}</text>`+original+plotted.map((s,i)=>{const x=20+(i%2)*width/2,y=height+18+Math.floor(i/2)*23;return `<line x1="${x}" x2="${x+18}" y1="${y-4}" y2="${y-4}" stroke="${s.color}" stroke-width="3" ${s.dash?`stroke-dasharray="${s.dash}"`:''}/><text x="${x+26}" y="${y}" fill="#183c38" font-family="Arial" font-size="10">${escapeHTML(s.label)}</text>`;}).join('')+`<text x="20" y="${height+37+rows*23}" fill="#697774" font-family="Arial" font-size="10">mortselect · Saved research point forecasts · Mortality rate × 100</text>`;downloadBlob(new XMLSerializer().serializeToString(copy),'image/svg+xml',`mortselect-${state.method}-${state.mode}.svg`);});
  $('share').addEventListener('click',async()=>{saveLocation();try{await navigator.clipboard.writeText(location.href);$('notice').textContent='View link copied. It includes your current filters.';}catch{$('notice').textContent='Your view is saved in the address bar. Copy that URL to share it.';}});
  $('chart').addEventListener('pointermove',event=>{if(!plotGeometry)return;const rect=$('chart').getBoundingClientRect();const pos=(event.clientX-rect.left)*plotGeometry.width/rect.width;const {pad,width,xs}=plotGeometry;inspect(Math.round((pos-pad.left)/(width-pad.left-pad.right)*(xs.length-1)));});
  $('chart').addEventListener('pointerleave',hideTooltip);$('chart').addEventListener('blur',hideTooltip);$('chart').addEventListener('focus',()=>inspect(focusIndex));
  $('chart').addEventListener('keydown',event=>{if(['ArrowLeft','ArrowRight','Home','End','Escape'].includes(event.key)){event.preventDefault();if(event.key==='Escape')return hideTooltip();inspect(event.key==='Home'?0:event.key==='End'?plotGeometry.xs.length-1:focusIndex+(event.key==='ArrowRight'?1:-1));}});
  new ResizeObserver(()=>{if(plotted.length){hideTooltip();drawChart();}}).observe($('chart-container'));
  window.addEventListener('hashchange',()=>{if(location.hash.startsWith('#view=')){readLocation();normalize();populateControls();loadMethod();}});
}
async function init(){
  setStatus('Loading the research dataset…');
  try{manifest=await fetchJSON('data/manifest.json');$('population-count').textContent=String(manifest.population_count).padStart(2,'0');$('method').innerHTML=manifest.methods.map(m=>`<option value="${m.id}">${escapeHTML(m.label)}</option>`).join('');$('method').disabled=false;readLocation();normalize();populateControls();await loadMethod();}
  catch{setStatus('The dataset could not be loaded. Check your connection and try again.',true);}
}
setupEvents();init();
