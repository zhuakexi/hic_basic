import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { LineMaterial } from "three/addons/lines/LineMaterial.js";
import { LineSegments2 } from "three/addons/lines/LineSegments2.js";
import { LineSegmentsGeometry } from "three/addons/lines/LineSegmentsGeometry.js";

const payloadNode = document.getElementById("scene-data");
const payload = JSON.parse(payloadNode.textContent);
payloadNode.remove();
const RENDERER_PROTOCOL_VERSION = 2;
const embeddedRenderConfig = payload.viewer.render_config || {
  schema_version: 1,
  view: {
    regions: payload.viewer.initial_regions || [],
    camera: { mode: "fit", fov: 34, up: [0, 1, 0] },
    section: { enabled: false, thickness: 2, offset: 0 },
  },
  layers: { show_particles: null, show_bonds: null, show_annotations: null },
  style: { profile: payload.viewer.export_profile || "molecular-studio-v1" },
  raster: { long_edge: 2400, width: null, height: null },
};

const viewport = document.getElementById("viewport");
const errorState = document.getElementById("error-state");
const particlesToggle = document.getElementById("particles-toggle");
const bondsToggle = document.getElementById("bonds-toggle");
const annotationsToggle = document.getElementById("annotations-toggle");
const pointSizeInput = document.getElementById("point-size");
const bondSizeInput = document.getElementById("bond-size");
const sectionToggle = document.getElementById("section-toggle");
const sectionThicknessInput = document.getElementById("section-thickness");
const sectionThicknessValue = document.getElementById("section-thickness-value");
const sectionOffsetInput = document.getElementById("section-offset");
const sectionOffsetValue = document.getElementById("section-offset-value");
const backgroundInput = document.getElementById("background-color");
const exportProfileSelect = document.getElementById("export-profile");
const fitButton = document.getElementById("fit-view");
const saveButton = document.getElementById("save-png");
const regionsButton = document.getElementById("regions-button");
const regionSummary = document.getElementById("region-summary");
const regionPanel = document.getElementById("region-panel");
const regionForm = document.getElementById("region-form");
const regionInput = document.getElementById("region-input");
const regionList = document.getElementById("region-list");
const regionError = document.getElementById("region-error");
const clearRegionsButton = document.getElementById("clear-regions");
const tooltip = document.getElementById("tooltip");
const selectionStatus = document.getElementById("selection-status");
const fpsStatus = document.getElementById("fps-status");
const rendererStatus = document.getElementById("renderer-status");

document.getElementById("dataset-name").textContent = payload.dataset.name;
document.getElementById("particle-count").textContent = `${payload.stats.particles.toLocaleString()} particles`;
document.getElementById("bond-count").textContent = `${payload.stats.bonds.toLocaleString()} bonds`;
document.getElementById("resolution").textContent = payload.dataset.resolution
  ? `${Math.round(payload.dataset.resolution / 1000)} kb`
  : "resolution unknown";

if (payload.viewer.background) backgroundInput.value = payload.viewer.background;
if (payload.viewer.particle_size) pointSizeInput.value = String(payload.viewer.particle_size);
if (payload.viewer.bond_size) bondSizeInput.value = String(payload.viewer.bond_size);
if (payload.viewer.export_profile) exportProfileSelect.value = payload.viewer.export_profile;
if (embeddedRenderConfig.layers.show_particles !== null) {
  particlesToggle.checked = embeddedRenderConfig.layers.show_particles;
}
if (embeddedRenderConfig.layers.show_bonds !== null) {
  bondsToggle.checked = embeddedRenderConfig.layers.show_bonds;
}
const embeddedSection = embeddedRenderConfig.view.section || {};
sectionToggle.checked = Boolean(embeddedSection.enabled);
if (embeddedSection.thickness !== undefined) {
  sectionThicknessInput.value = String(embeddedSection.thickness);
}
if (embeddedSection.offset !== undefined) sectionOffsetInput.value = String(embeddedSection.offset);

function decodeBase64(value, Type) {
  const raw = atob(value);
  const bytes = new Uint8Array(raw.length);
  const chunkSize = 1 << 18;
  for (let offset = 0; offset < raw.length; offset += chunkSize) {
    const end = Math.min(offset + chunkSize, raw.length);
    for (let index = offset; index < end; index += 1) {
      bytes[index] = raw.charCodeAt(index);
    }
  }
  return new Type(bytes.buffer);
}

const ARRAY_TYPES = {
  float32: Float32Array,
  uint8: Uint8Array,
  uint16: Uint16Array,
  uint32: Uint32Array,
};

function decodeSceneArray(value, legacyType) {
  if (typeof value === "string") return decodeBase64(value, legacyType);
  const Type = ARRAY_TYPES[value.dtype];
  if (!Type) throw new Error(`Unsupported scene array dtype: ${value.dtype}`);
  return decodeBase64(value.data, Type);
}

let renderer;
try {
  renderer = new THREE.WebGLRenderer({
    antialias: false,
    alpha: false,
    powerPreference: "default",
    preserveDrawingBuffer: true,
  });
} catch (error) {
  errorState.style.display = "grid";
  errorState.textContent = `WebGL 2 renderer could not start: ${error.message}`;
  throw error;
}

const gl = renderer.getContext();
const debugRendererInfo = gl.getExtension("WEBGL_debug_renderer_info");
const rendererLabel = debugRendererInfo
  ? gl.getParameter(debugRendererInfo.UNMASKED_RENDERER_WEBGL)
  : gl.getParameter(gl.RENDERER);
const softwareRenderer = /swiftshader|llvmpipe|software/i.test(rendererLabel);
renderer.setPixelRatio(softwareRenderer ? 0.72 : Math.min(window.devicePixelRatio || 1, 2));
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(backgroundInput.value, 1);
renderer.outputColorSpace = THREE.SRGBColorSpace;
renderer.localClippingEnabled = sectionToggle.checked;
viewport.appendChild(renderer.domElement);

const positions = decodeSceneArray(payload.arrays.positions, Float32Array);
const colors = decodeSceneArray(payload.arrays.colors, Uint8Array);
const chromosomeCodes = decodeSceneArray(payload.arrays.chromosome_codes, Uint8Array);
const genomicPositions = decodeSceneArray(payload.arrays.genomic_positions, Uint32Array);
const bondIndices = decodeSceneArray(payload.arrays.bond_indices, Uint32Array);
const annotationPayload = payload.annotations || { schema_version: 1, layers: [] };
if (annotationPayload.schema_version !== 1) {
  throw new Error(`Unsupported annotation schema version: ${annotationPayload.schema_version}`);
}
const annotationLayers = annotationPayload.layers.map((layer) => ({
  ...layer,
  annotations: layer.annotations.map((annotation) => ({
    ...annotation,
    target: {
      ...annotation.target,
      particleIndices: decodeSceneArray(annotation.target.particle_indices, Uint32Array),
      bondIndices: decodeSceneArray(annotation.target.bond_indices, Uint32Array),
    },
  })),
}));
const hasAnnotations = annotationLayers.some((layer) => layer.annotations.length > 0);
annotationsToggle.checked = embeddedRenderConfig.layers.show_annotations === null
  || embeddedRenderConfig.layers.show_annotations === undefined
  ? annotationLayers.some((layer) => layer.visible)
  : Boolean(embeddedRenderConfig.layers.show_annotations);
annotationsToggle.disabled = !hasAnnotations;

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(34, window.innerWidth / window.innerHeight, 0.01, 5000);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.075;
controls.rotateSpeed = 0.72;
controls.panSpeed = 0.7;
controls.zoomSpeed = 0.9;
controls.screenSpacePanning = true;

const positionAttribute = new THREE.BufferAttribute(positions, 3);
const colorAttribute = new THREE.BufferAttribute(colors, 3, true);
const monitorSectionPlanes = [new THREE.Plane(), new THREE.Plane()];

const particleGeometry = new THREE.BufferGeometry();
particleGeometry.setAttribute("position", positionAttribute);
particleGeometry.setAttribute("color", colorAttribute);
particleGeometry.computeBoundingBox();
particleGeometry.computeBoundingSphere();

const particleMaterial = new THREE.PointsMaterial({
  size: payload.viewer.base_point_size,
  sizeAttenuation: true,
  vertexColors: true,
  transparent: false,
  depthWrite: true,
  clippingPlanes: monitorSectionPlanes,
});
const particles = new THREE.Points(particleGeometry, particleMaterial);
particles.frustumCulled = true;
scene.add(particles);

const thinBondGeometry = new THREE.BufferGeometry();
thinBondGeometry.setAttribute("position", positionAttribute);
thinBondGeometry.setAttribute("color", colorAttribute);
thinBondGeometry.setIndex(new THREE.BufferAttribute(bondIndices, 1));
thinBondGeometry.computeBoundingBox();
thinBondGeometry.computeBoundingSphere();

const thinBondMaterial = new THREE.LineBasicMaterial({
  vertexColors: true,
  transparent: true,
  opacity: 0.58,
  depthWrite: true,
  clippingPlanes: monitorSectionPlanes,
});
const thinBonds = new THREE.LineSegments(thinBondGeometry, thinBondMaterial);
const bonds = new THREE.Group();
bonds.renderOrder = -1;
bonds.add(thinBonds);
scene.add(bonds);

const monitorAnnotations = new THREE.Group();
monitorAnnotations.renderOrder = 5;
scene.add(monitorAnnotations);
let monitorAnnotationResources = [];

let wideBonds = null;
let wideBondGeometry = null;
let wideBondMaterial = null;
let wideBondGeometryDirty = true;
let monitorBondIndices = bondIndices;

const haloGeometry = new THREE.SphereGeometry(payload.viewer.base_point_size * 1.9, 14, 10);
const haloMaterial = new THREE.MeshBasicMaterial({
  color: 0xd1495b,
  wireframe: true,
  transparent: true,
  opacity: 0.95,
  depthTest: false,
  clippingPlanes: monitorSectionPlanes,
});
const selectionHalo = new THREE.Mesh(haloGeometry, haloMaterial);
selectionHalo.visible = false;
selectionHalo.renderOrder = 4;
scene.add(selectionHalo);

const allSelection = {
  name: "All genome",
  particleCount: payload.stats.particles,
  bondCount: payload.stats.bonds,
  bbox: payload.stats.bbox,
  particleIndices: null,
  bondIndices,
};

const groupByName = new Map(payload.groups.map((group) => [group.name.toLowerCase(), group]));
const groupByCode = new Map(payload.groups.map((group) => [group.code, group]));
let activeRegions = [];
let activeSelection = allSelection;
let controlsActive = false;
let lastPickTime = 0;
let lastFrameTime = performance.now();
let fpsWindowStart = lastFrameTime;
let fpsFrames = 0;
let firstFrameTime = null;

function sectionState(sectionConfig, targetCamera, target) {
  const enabled = Boolean(sectionConfig && sectionConfig.enabled);
  const thickness = finiteNumber(sectionConfig && sectionConfig.thickness === undefined
    ? 2
    : sectionConfig.thickness, "section thickness");
  const offset = finiteNumber(sectionConfig && sectionConfig.offset === undefined
    ? 0
    : sectionConfig.offset, "section offset");
  if (!(thickness > 0)) throw new Error("section thickness must be greater than 0");
  const normal = targetCamera.position.clone().sub(target).normalize();
  const center = target.clone().addScaledVector(normal, offset);
  const half = thickness * 0.5;
  const lower = center.clone().addScaledVector(normal, -half);
  const upper = center.clone().addScaledVector(normal, half);
  return {
    enabled,
    thickness,
    offset,
    normal,
    center,
    planes: [
      new THREE.Plane().setFromNormalAndCoplanarPoint(normal, lower),
      new THREE.Plane().setFromNormalAndCoplanarPoint(normal.clone().negate(), upper),
    ],
  };
}

function updateMonitorSectionPlanes() {
  const state = sectionState(
    {
      enabled: sectionToggle.checked,
      thickness: sectionThicknessInput.value,
      offset: sectionOffsetInput.value,
    },
    camera,
    controls.target,
  );
  monitorSectionPlanes[0].copy(state.planes[0]);
  monitorSectionPlanes[1].copy(state.planes[1]);
  renderer.localClippingEnabled = state.enabled;
  return state;
}

function configureSectionInputs(selection = activeSelection) {
  const currentThickness = Number(sectionThicknessInput.value);
  const currentOffset = Number(sectionOffsetInput.value);
  const minimum = new THREE.Vector3(...selection.bbox[0]);
  const maximum = new THREE.Vector3(...selection.bbox[1]);
  const diagonal = Math.max(maximum.distanceTo(minimum), 1e-3);
  const minimumThickness = Math.max(diagonal / 500, 1e-4);
  const maximumThickness = Math.max(diagonal * 1.25, Number(embeddedSection.thickness) || 0);
  sectionThicknessInput.min = String(minimumThickness);
  sectionThicknessInput.max = String(maximumThickness);
  sectionThicknessInput.step = "any";
  sectionOffsetInput.min = String(-diagonal);
  sectionOffsetInput.max = String(diagonal);
  sectionOffsetInput.step = "any";
  sectionThicknessInput.value = String(THREE.MathUtils.clamp(
    currentThickness,
    minimumThickness,
    maximumThickness,
  ));
  sectionOffsetInput.value = String(THREE.MathUtils.clamp(currentOffset, -diagonal, diagonal));
}

function formatSectionValue(value) {
  const numeric = Number(value);
  if (Math.abs(numeric) >= 100) return numeric.toFixed(0);
  if (Math.abs(numeric) >= 10) return numeric.toFixed(1);
  return numeric.toFixed(2).replace(/\.00$/, "");
}

function updateSectionControls() {
  document.querySelectorAll(".field-section").forEach((field) => {
    field.dataset.disabled = String(!sectionToggle.checked);
  });
  sectionThicknessInput.disabled = !sectionToggle.checked;
  sectionOffsetInput.disabled = !sectionToggle.checked;
  sectionThicknessValue.textContent = formatSectionValue(sectionThicknessInput.value);
  sectionOffsetValue.textContent = formatSectionValue(sectionOffsetInput.value);
  updateMonitorSectionPlanes();
}

function lowerBound(array, value, start, end) {
  let left = start;
  let right = end;
  while (left < right) {
    const middle = (left + right) >>> 1;
    if (array[middle] < value) left = middle + 1;
    else right = middle;
  }
  return left;
}

function upperBound(array, value, start, end) {
  let left = start;
  let right = end;
  while (left < right) {
    const middle = (left + right) >>> 1;
    if (array[middle] <= value) left = middle + 1;
    else right = middle;
  }
  return left;
}

function parseCoordinate(value) {
  const match = value.trim().match(/^([0-9][0-9,_.]*(?:\.[0-9]+)?)\s*(bp|kb|mb|gb)?$/i);
  if (!match) throw new Error(`Invalid coordinate: ${value}`);
  const numeric = Number(match[1].replace(/[,_]/g, ""));
  const multiplier = { bp: 1, kb: 1e3, mb: 1e6, gb: 1e9 }[(match[2] || "bp").toLowerCase()];
  const result = Math.round(numeric * multiplier);
  if (!Number.isSafeInteger(result) || result < 0) throw new Error(`Coordinate out of range: ${value}`);
  return result;
}

function resolveChromosomes(rawName) {
  const stripped = rawName.trim();
  const alternatives = [stripped];
  if (!/^chr/i.test(stripped)) alternatives.push(`chr${stripped}`);

  for (const alternative of alternatives) {
    const exact = groupByName.get(alternative.toLowerCase());
    if (exact) return [exact];
  }

  for (const alternative of alternatives) {
    const prefix = `${alternative.toLowerCase()}(`;
    const homologs = payload.groups.filter((group) => group.name.toLowerCase().startsWith(prefix));
    if (homologs.length) return homologs;
  }
  throw new Error(`Chromosome not found: ${rawName}`);
}

function parseRegion(rawRegion) {
  const value = rawRegion.trim();
  if (!value) throw new Error("Enter a genome region");
  const colon = value.indexOf(":");
  const chromosomeToken = colon === -1 ? value : value.slice(0, colon).trim();
  const coordinateToken = colon === -1 ? "" : value.slice(colon + 1).trim();
  const groups = resolveChromosomes(chromosomeToken);
  let start = 0;
  let end = Number.MAX_SAFE_INTEGER;
  if (coordinateToken) {
    const rangeMatch = coordinateToken.match(/^(.+?)\s*-\s*(.+)$/);
    if (!rangeMatch) throw new Error("Use chr:start-end for a partial region");
    start = parseCoordinate(rangeMatch[1]);
    end = parseCoordinate(rangeMatch[2]);
    if (start > end) throw new Error("Region start is greater than its end");
  }
  return { raw: value, groups, start, end };
}

function mergeIntervals(regions) {
  const byCode = new Map();
  regions.forEach((region) => {
    region.groups.forEach((group) => {
      const intervals = byCode.get(group.code) || [];
      intervals.push([region.start, region.end]);
      byCode.set(group.code, intervals);
    });
  });

  const merged = new Map();
  byCode.forEach((intervals, code) => {
    intervals.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    const result = [];
    intervals.forEach(([start, end]) => {
      const previous = result[result.length - 1];
      if (previous && start <= previous[1] + 1) previous[1] = Math.max(previous[1], end);
      else result.push([start, end]);
    });
    merged.set(code, result);
  });
  return merged;
}

function selectionFromRegions(regions) {
  if (!regions.length) return allSelection;
  const merged = mergeIntervals(regions);
  const spans = [];
  Array.from(merged.entries()).sort((a, b) => a[0] - b[0]).forEach(([code, intervals]) => {
    const group = groupByCode.get(code);
    const groupEnd = group.start + group.count;
    intervals.forEach(([regionStart, regionEnd]) => {
      const start = lowerBound(genomicPositions, regionStart, group.start, groupEnd);
      const end = upperBound(genomicPositions, regionEnd, group.start, groupEnd);
      if (end > start) spans.push([start, end]);
    });
  });
  const particleCount = spans.reduce((total, [start, end]) => total + end - start, 0);
  if (!particleCount) throw new Error("No particles overlap the selected regions");

  const particleIndices = new Uint32Array(particleCount);
  let particleOffset = 0;
  spans.forEach(([start, end]) => {
    for (let index = start; index < end; index += 1) particleIndices[particleOffset++] = index;
  });

  const membership = new Uint8Array(payload.stats.particles);
  for (let index = 0; index < particleCount; index += 1) membership[particleIndices[index]] = 1;
  const selectedBonds = [];
  for (let index = 0; index < bondIndices.length; index += 2) {
    const left = bondIndices[index];
    const right = bondIndices[index + 1];
    if (membership[left] && membership[right]) selectedBonds.push(left, right);
  }

  const selectedBondIndices = Uint32Array.from(selectedBonds);
  const minimum = [Infinity, Infinity, Infinity];
  const maximum = [-Infinity, -Infinity, -Infinity];
  for (let offset = 0; offset < particleIndices.length; offset += 1) {
    const particleIndex = particleIndices[offset] * 3;
    for (let axis = 0; axis < 3; axis += 1) {
      const value = positions[particleIndex + axis];
      minimum[axis] = Math.min(minimum[axis], value);
      maximum[axis] = Math.max(maximum[axis], value);
    }
  }
  return {
    name: regions.length === 1 ? regions[0].raw : `${regions.length} regions`,
    particleCount,
    bondCount: selectedBondIndices.length / 2,
    bbox: [minimum, maximum],
    particleIndices,
    bondIndices: selectedBondIndices,
  };
}

function profileColorArray(profileName, targetCamera, selection) {
  const profile = EXPORT_PROFILES[profileName];
  const output = new Float32Array(colors.length);
  const color = new THREE.Color();
  const hsl = {};

  for (let particleIndex = 0; particleIndex < chromosomeCodes.length; particleIndex += 1) {
    const index = particleIndex * 3;
    color.setRGB(colors[index] / 255, colors[index + 1] / 255, colors[index + 2] / 255);
    color.getHSL(hsl);
    hsl.s = THREE.MathUtils.clamp(
      hsl.s * profile.saturationScale,
      profile.saturationRange[0],
      profile.saturationRange[1],
    );
    hsl.l = THREE.MathUtils.clamp(
      hsl.l * profile.lightnessScale,
      profile.lightnessRange[0],
      profile.lightnessRange[1],
    );
    color.setHSL(hsl.h, hsl.s, hsl.l);
    output[index] = color.r;
    output[index + 1] = color.g;
    output[index + 2] = color.b;
  }

  targetCamera.updateMatrixWorld(true);
  targetCamera.matrixWorldInverse.copy(targetCamera.matrixWorld).invert();
  const matrix = targetCamera.matrixWorldInverse.elements;
  const selected = selection.particleIndices;
  const count = selected ? selected.length : payload.stats.particles;
  let nearDepth = Infinity;
  let farDepth = -Infinity;
  for (let offset = 0; offset < count; offset += 1) {
    const particle = selected ? selected[offset] : offset;
    const index = particle * 3;
    const depth = -(matrix[2] * positions[index]
      + matrix[6] * positions[index + 1]
      + matrix[10] * positions[index + 2]
      + matrix[14]);
    nearDepth = Math.min(nearDepth, depth);
    farDepth = Math.max(farDepth, depth);
  }
  const depthSpan = Math.max(farDepth - nearDepth, 1e-6);
  for (let offset = 0; offset < count; offset += 1) {
    const particle = selected ? selected[offset] : offset;
    const index = particle * 3;
    const depth = -(matrix[2] * positions[index]
      + matrix[6] * positions[index + 1]
      + matrix[10] * positions[index + 2]
      + matrix[14]);
    const depthFraction = THREE.MathUtils.clamp((depth - nearDepth) / depthSpan, 0, 1);
    const shade = 1 - profile.depthCue * depthFraction;
    output[index] *= shade;
    output[index + 1] *= shade;
    output[index + 2] *= shade;
  }
  return output;
}

function annotationColorWriter(mark, profile, targetCamera, selection) {
  const color = new THREE.Color(mark.color);
  if (mark.material === "unlit") {
    return (particleIndex, target, offset) => {
      target[offset] = color.r;
      target[offset + 1] = color.g;
      target[offset + 2] = color.b;
    };
  }

  const hsl = {};
  color.getHSL(hsl);
  hsl.s = THREE.MathUtils.clamp(
    hsl.s * profile.saturationScale,
    profile.saturationRange[0],
    profile.saturationRange[1],
  );
  hsl.l = THREE.MathUtils.clamp(
    hsl.l * profile.lightnessScale,
    profile.lightnessRange[0],
    profile.lightnessRange[1],
  );
  color.setHSL(hsl.h, hsl.s, hsl.l);
  targetCamera.updateMatrixWorld(true);
  targetCamera.matrixWorldInverse.copy(targetCamera.matrixWorld).invert();
  const matrix = targetCamera.matrixWorldInverse.elements;
  const selected = selection.particleIndices;
  const count = selected ? selected.length : payload.stats.particles;
  let nearDepth = Infinity;
  let farDepth = -Infinity;
  const particleDepth = (particle) => {
    const index = particle * 3;
    return -(matrix[2] * positions[index]
      + matrix[6] * positions[index + 1]
      + matrix[10] * positions[index + 2]
      + matrix[14]);
  };
  for (let offset = 0; offset < count; offset += 1) {
    const depth = particleDepth(selected ? selected[offset] : offset);
    nearDepth = Math.min(nearDepth, depth);
    farDepth = Math.max(farDepth, depth);
  }
  const depthSpan = Math.max(farDepth - nearDepth, 1e-6);
  return (particleIndex, target, offset) => {
    const depthFraction = THREE.MathUtils.clamp(
      (particleDepth(particleIndex) - nearDepth) / depthSpan,
      0,
      1,
    );
    const shade = 1 - profile.depthCue * depthFraction;
    target[offset] = color.r * shade;
    target[offset + 1] = color.g * shade;
    target[offset + 2] = color.b * shade;
  };
}

function lineArrays(indices, colorSource = colors, normalizedBytes = true) {
  const linePositions = new Float32Array(indices.length * 3);
  const lineColors = new Float32Array(indices.length * 3);
  for (let endpoint = 0; endpoint < indices.length; endpoint += 1) {
    const particleIndex = indices[endpoint] * 3;
    const targetIndex = endpoint * 3;
    linePositions[targetIndex] = positions[particleIndex];
    linePositions[targetIndex + 1] = positions[particleIndex + 1];
    linePositions[targetIndex + 2] = positions[particleIndex + 2];
    const scale = normalizedBytes ? 1 / 255 : 1;
    lineColors[targetIndex] = colorSource[particleIndex] * scale;
    lineColors[targetIndex + 1] = colorSource[particleIndex + 1] * scale;
    lineColors[targetIndex + 2] = colorSource[particleIndex + 2] * scale;
  }
  return { linePositions, lineColors };
}

function annotationBondEndpoints(annotation, selection) {
  let membership = null;
  if (selection.particleIndices) {
    membership = new Uint8Array(payload.stats.particles);
    for (let index = 0; index < selection.particleIndices.length; index += 1) {
      membership[selection.particleIndices[index]] = 1;
    }
  }
  const endpoints = [];
  for (let index = 0; index < annotation.target.bondIndices.length; index += 1) {
    const bondIndex = annotation.target.bondIndices[index] * 2;
    const left = bondIndices[bondIndex];
    const right = bondIndices[bondIndex + 1];
    if (!membership || (membership[left] && membership[right])) endpoints.push(left, right);
  }
  return Uint32Array.from(endpoints);
}

function baseMonitorBondEndpoints(selection) {
  if (!annotationsToggle.checked || !hasAnnotations) return selection.bondIndices;
  const annotated = new Set();
  annotationLayers.forEach((layer) => {
    if (!layer.visible) return;
    layer.annotations.forEach((annotation) => {
      const endpoints = annotationBondEndpoints(annotation, selection);
      for (let index = 0; index < endpoints.length; index += 2) {
        annotated.add(`${endpoints[index]}:${endpoints[index + 1]}`);
      }
    });
  });
  if (!annotated.size) return selection.bondIndices;
  const remaining = [];
  for (let index = 0; index < selection.bondIndices.length; index += 2) {
    const left = selection.bondIndices[index];
    const right = selection.bondIndices[index + 1];
    if (!annotated.has(`${left}:${right}`)) remaining.push(left, right);
  }
  return Uint32Array.from(remaining);
}

function updateMonitorBondSelection() {
  monitorBondIndices = baseMonitorBondEndpoints(activeSelection);
  thinBondGeometry.setIndex(new THREE.BufferAttribute(monitorBondIndices, 1));
  wideBondGeometryDirty = true;
  applyBondSize();
}

function disposeMonitorAnnotations() {
  monitorAnnotations.clear();
  monitorAnnotationResources.forEach(({ geometry, materials }) => {
    geometry.dispose();
    materials.forEach((material) => material.dispose());
  });
  monitorAnnotationResources = [];
}

function rebuildMonitorAnnotations() {
  disposeMonitorAnnotations();
  monitorAnnotations.visible = annotationsToggle.checked && hasAnnotations;
  if (!monitorAnnotations.visible) return;
  const resolution = new THREE.Vector2(window.innerWidth, window.innerHeight);
  const scale = Math.max(window.innerWidth / 1400, 0.5);
  annotationLayers.forEach((layer, layerIndex) => {
    if (!layer.visible) return;
    layer.annotations.forEach((annotation, annotationIndex) => {
      const endpoints = annotationBondEndpoints(annotation, activeSelection);
      if (!endpoints.length) return;
      annotation.marks.forEach((mark, markIndex) => {
        if (mark.kind !== "tube") throw new Error(`Unsupported annotation mark: ${mark.kind}`);
        const arrays = lineArrays(endpoints);
        const geometry = new LineSegmentsGeometry();
        geometry.setPositions(arrays.linePositions);
        const depthTest = mark.occlusion === "depth_tested";
        const renderOrder = 10 + layerIndex * 100 + annotationIndex * 4 + markIndex * 2;
        const materials = [];
        if (mark.outline_width > mark.tube_width && mark.outline_color) {
          const outline = makeLineLayer(geometry, mark.outline_width * scale, resolution, {
            vertexColors: false,
            color: mark.outline_color,
            opacity: mark.opacity,
            transparent: true,
            depthTest,
            depthWrite: false,
            toneMapped: mark.material !== "unlit",
            renderOrder,
            clippingPlanes: monitorSectionPlanes,
          });
          monitorAnnotations.add(outline.line);
          materials.push(outline.material);
        }
        const tube = makeLineLayer(geometry, mark.tube_width * scale, resolution, {
          vertexColors: false,
          color: mark.color,
          opacity: mark.opacity,
          transparent: true,
          depthTest,
          depthWrite: false,
          toneMapped: mark.material !== "unlit",
          renderOrder: renderOrder + 1,
          clippingPlanes: monitorSectionPlanes,
        });
        monitorAnnotations.add(tube.line);
        materials.push(tube.material);
        monitorAnnotationResources.push({ geometry, materials });
      });
    });
  });
}

function disposeWideBonds() {
  if (wideBonds) bonds.remove(wideBonds);
  if (wideBondGeometry) wideBondGeometry.dispose();
  if (wideBondMaterial) wideBondMaterial.dispose();
  wideBonds = null;
  wideBondGeometry = null;
  wideBondMaterial = null;
  wideBondGeometryDirty = true;
}

function rebuildWideBonds() {
  disposeWideBonds();
  const arrays = lineArrays(monitorBondIndices);
  wideBondGeometry = new LineSegmentsGeometry();
  wideBondGeometry.setPositions(arrays.linePositions);
  wideBondGeometry.setColors(arrays.lineColors);
  wideBondMaterial = new LineMaterial({
    linewidth: Number(bondSizeInput.value),
    vertexColors: true,
    transparent: true,
    opacity: 0.70,
    depthWrite: true,
    alphaToCoverage: true,
    clippingPlanes: monitorSectionPlanes,
  });
  wideBondMaterial.resolution.set(window.innerWidth, window.innerHeight);
  wideBonds = new LineSegments2(wideBondGeometry, wideBondMaterial);
  wideBonds.renderOrder = -1;
  bonds.add(wideBonds);
  wideBondGeometryDirty = false;
}

function applyBondSize() {
  const size = Number(bondSizeInput.value);
  if (size <= 1.05) {
    thinBonds.visible = true;
    disposeWideBonds();
    return;
  }
  if (wideBondGeometryDirty || !wideBonds) rebuildWideBonds();
  wideBondMaterial.linewidth = size;
  wideBonds.visible = true;
  thinBonds.visible = false;
}

function applySelection(selection, fit = true) {
  activeSelection = selection;
  if (selection.particleIndices) {
    particleGeometry.setIndex(new THREE.BufferAttribute(selection.particleIndices, 1));
  } else {
    particleGeometry.setIndex(null);
  }
  selectionHalo.visible = false;
  selectionStatus.textContent = `${selection.name}: ${selection.particleCount.toLocaleString()} particles, ${selection.bondCount.toLocaleString()} bonds`;
  updateMonitorBondSelection();
  rebuildMonitorAnnotations();
  if (fit) fitView(selection);
  configureSectionInputs(selection);
  updateSectionControls();
}

function frameCamera(selection, targetCamera, target, directionOverride = null) {
  const [minimum, maximum] = selection.bbox;
  const min = new THREE.Vector3(...minimum);
  const max = new THREE.Vector3(...maximum);
  const box = new THREE.Box3(min, max);
  const center = box.getCenter(new THREE.Vector3());
  const size = box.getSize(new THREE.Vector3());
  const radius = Math.max(size.length() * 0.5, 1);
  const currentDirection = targetCamera.position.clone().sub(target);
  const direction = directionOverride || (currentDirection.lengthSq() > 0
    ? currentDirection.normalize()
    : new THREE.Vector3(1.08, 0.76, 1.32).normalize());
  const verticalFov = THREE.MathUtils.degToRad(targetCamera.fov);
  const horizontalFov = 2 * Math.atan(Math.tan(verticalFov / 2) * targetCamera.aspect);
  const limitingFov = Math.min(verticalFov, horizontalFov);
  const distance = (radius / Math.sin(limitingFov / 2)) * 0.78;

  target.copy(center);
  targetCamera.position.copy(center).addScaledVector(direction, distance);
  targetCamera.near = Math.max(radius / 250, 0.01);
  targetCamera.far = radius * 20;
  targetCamera.lookAt(target);
  targetCamera.updateProjectionMatrix();
}

function fitView(selection = activeSelection, targetCamera = camera, targetControls = controls) {
  frameCamera(selection, targetCamera, targetControls.target);
  targetControls.update();
}

function renderRegionList() {
  regionList.replaceChildren();
  activeRegions.forEach((region, index) => {
    const chip = document.createElement("span");
    chip.className = "region-chip";
    const label = document.createElement("span");
    label.textContent = region.raw;
    const remove = document.createElement("button");
    remove.type = "button";
    remove.textContent = "\u00d7";
    remove.title = `Remove ${region.raw}`;
    remove.setAttribute("aria-label", `Remove ${region.raw}`);
    remove.addEventListener("click", () => {
      activeRegions.splice(index, 1);
      applyRegions();
    });
    chip.append(label, remove);
    regionList.appendChild(chip);
  });
  regionSummary.textContent = activeRegions.length ? String(activeRegions.length) : "All";
  clearRegionsButton.disabled = activeRegions.length === 0;
}

function applyRegions() {
  try {
    const selection = selectionFromRegions(activeRegions);
    regionError.textContent = "";
    renderRegionList();
    applySelection(selection);
  } catch (error) {
    regionError.textContent = error.message;
  }
}

function addRegions(rawInput) {
  const values = rawInput.split(/[;\n]+/).map((value) => value.trim()).filter(Boolean);
  if (!values.length) throw new Error("Enter a genome region");
  const parsed = values.map(parseRegion);
  const candidate = activeRegions.concat(parsed);
  selectionFromRegions(candidate);
  activeRegions = candidate;
  regionInput.value = "";
  applyRegions();
}

regionsButton.addEventListener("click", () => {
  const opening = regionPanel.hidden;
  regionPanel.hidden = !opening;
  regionsButton.setAttribute("aria-expanded", String(opening));
  if (opening) {
    positionRegionPanel();
    regionInput.focus();
  }
});

regionForm.addEventListener("submit", (event) => {
  event.preventDefault();
  try {
    addRegions(regionInput.value);
    regionError.textContent = "";
  } catch (error) {
    regionError.textContent = error.message;
  }
});

clearRegionsButton.addEventListener("click", () => {
  activeRegions = [];
  applyRegions();
});

document.addEventListener("pointerdown", (event) => {
  if (!regionPanel.hidden && !regionPanel.contains(event.target) && !regionsButton.contains(event.target)) {
    regionPanel.hidden = true;
    regionsButton.setAttribute("aria-expanded", "false");
  }
});

particlesToggle.addEventListener("change", () => {
  particles.visible = particlesToggle.checked;
});

bondsToggle.addEventListener("change", () => {
  bonds.visible = bondsToggle.checked;
});

annotationsToggle.addEventListener("change", () => {
  updateMonitorBondSelection();
  rebuildMonitorAnnotations();
});

pointSizeInput.addEventListener("input", () => {
  particleMaterial.size = payload.viewer.base_point_size * Number(pointSizeInput.value);
  selectionHalo.scale.setScalar(Number(pointSizeInput.value));
});

bondSizeInput.addEventListener("input", applyBondSize);
sectionToggle.addEventListener("change", updateSectionControls);
sectionThicknessInput.addEventListener("input", updateSectionControls);
sectionOffsetInput.addEventListener("input", updateSectionControls);

backgroundInput.addEventListener("input", () => {
  renderer.setClearColor(backgroundInput.value, 1);
});

fitButton.addEventListener("click", () => fitView(activeSelection));

const EXPORT_PROFILES = {
  "section-pymol-v1": {
    label: "section-pymol-v1",
    mode: "layered-section",
    background: "#ffffff",
    innerWidth: 4.6,
    outlineWidth: 4.6,
    highlightWidth: 0.1,
    pointSize: 2.2,
    outline: 0xffffff,
    outlineOpacity: 0,
    highlightOpacity: 0,
    saturationScale: 1.0,
    saturationRange: [0, 1],
    lightnessScale: 1.0,
    lightnessRange: [0, 1],
    depthCue: 0,
    exposure: 1.0,
    defaultParticles: false,
    showColorbar: true,
  },
  "section-studio-v1": {
    label: "section-studio-v1",
    mode: "layered-section",
    background: "#f7f9f7",
    innerWidth: 4.0,
    outlineWidth: 4.0,
    highlightWidth: 0.72,
    pointSize: 2.5,
    outline: 0x263335,
    outlineOpacity: 0,
    highlightOpacity: 0,
    saturationScale: 0.92,
    saturationRange: [0.18, 0.92],
    lightnessScale: 0.93,
    lightnessRange: [0.20, 0.78],
    depthCue: 0.08,
    exposure: 1.02,
    defaultParticles: false,
    showColorbar: true,
  },
  "xu2025-v1": {
    label: "xu2025-v1",
    mode: "layered-legacy",
    background: "#ffffff",
    innerWidth: 3.0,
    outlineWidth: 4.0,
    highlightWidth: 0.76,
    pointSize: 2.4,
    outline: 0x303638,
    outlineOpacity: 0.32,
    highlightOpacity: 0.34,
    saturationScale: 0.62,
    saturationRange: [0, 1],
    lightnessScale: 0.72,
    lightnessRange: [0.28, 0.57],
    depthCue: 0.30,
    exposure: 0.98,
  },
  "publication-v1": {
    label: "publication-v1",
    mode: "layered-legacy",
    background: "#f5f6f3",
    innerWidth: 3.4,
    outlineWidth: 4.6,
    highlightWidth: 1.02,
    pointSize: 2.8,
    outline: 0x20292c,
    outlineOpacity: 0.38,
    highlightOpacity: 0.42,
    saturationScale: 0.82,
    saturationRange: [0.42, 0.78],
    lightnessScale: 0.76,
    lightnessRange: [0.30, 0.58],
    depthCue: 0.38,
    exposure: 1.04,
  },
  "xu2025-v2": {
    label: "xu2025-v2",
    mode: "mesh-tube",
    background: "#ffffff",
    tubeWidth: 3.35,
    silhouetteWidth: 4.0,
    silhouette: 0x4a5152,
    radialSegments: 5,
    saturationScale: 0.60,
    saturationRange: [0.34, 0.68],
    lightnessScale: 0.78,
    lightnessRange: [0.32, 0.61],
    depthCue: 0.22,
    exposure: 1.02,
    lighting: {
      hemisphere: 1.35,
      key: 2.40,
      fill: 0.72,
      rim: 1.05,
      roughness: 0.68,
    },
  },
  "molecular-studio-v1": {
    label: "molecular-studio-v1",
    mode: "mesh-tube",
    background: "#111716",
    tubeWidth: 3.65,
    silhouetteWidth: 4.15,
    silhouette: 0x050807,
    radialSegments: 5,
    saturationScale: 0.86,
    saturationRange: [0.48, 0.82],
    lightnessScale: 0.82,
    lightnessRange: [0.36, 0.64],
    depthCue: 0.18,
    exposure: 1.08,
    lighting: {
      hemisphere: 0.82,
      key: 3.15,
      fill: 0.95,
      rim: 1.70,
      roughness: 0.48,
    },
  },
};

function cloneJson(value) {
  return JSON.parse(JSON.stringify(value));
}

function finiteNumber(value, name) {
  const result = Number(value);
  if (!Number.isFinite(result)) throw new Error(`${name} must be a finite number`);
  return result;
}

function positiveWidth(value, fallback, name) {
  const result = value === null || value === undefined ? fallback : finiteNumber(value, name);
  if (!(result > 0 && result <= 64)) throw new Error(`${name} must be greater than 0 and at most 64`);
  return result;
}

function normalizeRenderConfig(config) {
  const candidate = cloneJson(config || embeddedRenderConfig);
  if (candidate.schema_version !== 1) {
    throw new Error(`Unsupported render config schema version: ${candidate.schema_version}`);
  }
  const view = candidate.view || {};
  const layers = candidate.layers || {};
  const style = candidate.style || {};
  const raster = candidate.raster || {};
  const profileName = style.profile || exportProfileSelect.value;
  const profile = EXPORT_PROFILES[profileName];
  if (!profile) throw new Error(`Unknown export profile: ${profileName}`);

  const cameraConfig = view.camera || { mode: "fit" };
  const cameraMode = cameraConfig.mode || "fit";
  if (cameraMode !== "fit" && cameraMode !== "explicit") {
    throw new Error(`Unknown camera mode: ${cameraMode}`);
  }
  const cameraState = {
    mode: cameraMode,
    position: cameraConfig.position || null,
    target: cameraConfig.target || null,
    up: cameraConfig.up || [0, 1, 0],
    fov: finiteNumber(cameraConfig.fov === undefined ? 34 : cameraConfig.fov, "camera fov"),
    direction: cameraConfig.direction || null,
  };
  if (cameraState.direction) {
    if (!Array.isArray(cameraState.direction) || cameraState.direction.length !== 3) {
      throw new Error("Camera direction must contain three numbers");
    }
    cameraState.direction = cameraState.direction.map((value) => finiteNumber(value, "camera direction"));
    if (new THREE.Vector3(...cameraState.direction).lengthSq() <= 1e-20) {
      throw new Error("Camera direction must be non-zero");
    }
  }
  if (cameraMode === "explicit" && (!cameraState.position || !cameraState.target)) {
    throw new Error("Explicit camera requires position and target");
  }
  if (cameraMode === "explicit" && cameraState.direction) {
    throw new Error("Explicit camera does not accept direction");
  }
  const sectionConfig = view.section || {};
  const sectionThickness = finiteNumber(
    sectionConfig.thickness === undefined ? 2 : sectionConfig.thickness,
    "section thickness",
  );
  if (!(sectionThickness > 0)) throw new Error("section thickness must be greater than 0");
  const section = {
    enabled: Boolean(sectionConfig.enabled),
    thickness: sectionThickness,
    offset: finiteNumber(sectionConfig.offset === undefined ? 0 : sectionConfig.offset, "section offset"),
  };

  const width = raster.width === null || raster.width === undefined ? null : Number(raster.width);
  const height = raster.height === null || raster.height === undefined ? null : Number(raster.height);
  if ((width === null) !== (height === null)) throw new Error("Raster width and height must be provided together");
  const longEdge = Number(raster.long_edge || raster.longEdge || 2400);
  const profileParticleDefault = profile.pointSize || 3.0;
  const particlesSpecified = layers.show_particles !== null && layers.show_particles !== undefined;
  const showParticles = !particlesSpecified
    ? Boolean(profile.defaultParticles === undefined
      ? profile.mode === "layered-legacy"
      : profile.defaultParticles)
    : Boolean(layers.show_particles);
  const showBonds = layers.show_bonds === null || layers.show_bonds === undefined
    ? true
    : Boolean(layers.show_bonds);
  const showAnnotations = layers.show_annotations === null || layers.show_annotations === undefined
    ? annotationLayers.some((layer) => layer.visible)
    : Boolean(layers.show_annotations);

  return {
    schema_version: 1,
    view: {
      regions: Array.isArray(view.regions) ? view.regions.map(String) : [],
      camera: cameraState,
      section,
    },
    layers: {
      show_particles: showParticles,
      show_bonds: showBonds,
      show_annotations: showAnnotations,
    },
    particlesSpecified,
    style: {
      profile: profileName,
      background: style.background || profile.background,
      tube_width: positiveWidth(
        style.tube_width,
        profile.mode.startsWith("layered") ? profile.innerWidth : profile.tubeWidth,
        "tube_width",
      ),
      silhouette_width: positiveWidth(
        style.silhouette_width,
        profile.mode.startsWith("layered") ? profile.outlineWidth : profile.silhouetteWidth,
        "silhouette_width",
      ),
      particle_size: positiveWidth(style.particle_size, profileParticleDefault, "particle_size"),
      show_colorbar: style.show_colorbar === null || style.show_colorbar === undefined
        ? Boolean(profile.showColorbar)
        : Boolean(style.show_colorbar),
    },
    raster: { long_edge: longEdge, width, height },
    profile,
  };
}

function interactiveRenderConfig() {
  const config = cloneJson(embeddedRenderConfig);
  config.view.regions = activeRegions.map((region) => region.raw);
  config.view.camera = {
    mode: "explicit",
    position: camera.position.toArray(),
    target: controls.target.toArray(),
    up: camera.up.toArray(),
    fov: camera.fov,
  };
  config.view.section = {
    enabled: sectionToggle.checked,
    thickness: Number(sectionThicknessInput.value),
    offset: Number(sectionOffsetInput.value),
  };
  config.layers = {
    show_particles: particlesToggle.checked,
    show_bonds: bondsToggle.checked,
    show_annotations: annotationsToggle.checked,
  };
  config.style.profile = exportProfileSelect.value;
  return config;
}

function legacyRenderConfig(profileName, options) {
  const config = interactiveRenderConfig();
  config.style.profile = profileName || exportProfileSelect.value;
  if (options.longEdge !== undefined) config.raster.long_edge = options.longEdge;
  if (options.width !== undefined) config.raster.width = options.width;
  if (options.height !== undefined) config.raster.height = options.height;
  return config;
}

function selectionForRegions(regionStrings) {
  const parsed = regionStrings.map(parseRegion);
  return { regions: parsed, selection: selectionFromRegions(parsed) };
}

function configuredCamera(cameraConfig, selection, aspect) {
  const targetCamera = new THREE.PerspectiveCamera(cameraConfig.fov, aspect, 0.01, 5000);
  targetCamera.up.fromArray(cameraConfig.up);
  const target = new THREE.Vector3();
  if (cameraConfig.mode === "explicit") {
    targetCamera.position.fromArray(cameraConfig.position);
    target.fromArray(cameraConfig.target);
    const extent = new THREE.Vector3(...selection.bbox[1]).sub(new THREE.Vector3(...selection.bbox[0]));
    const radius = Math.max(extent.length() * 0.5, 1);
    targetCamera.near = Math.max(radius / 250, 0.01);
    targetCamera.far = Math.max(radius * 20, targetCamera.position.distanceTo(target) + radius * 4);
    targetCamera.lookAt(target);
    targetCamera.updateProjectionMatrix();
  } else {
    const fitDirection = cameraConfig.direction
      ? new THREE.Vector3(...cameraConfig.direction).normalize()
      : new THREE.Vector3(1.08, 0.76, 1.32).normalize();
    frameCamera(
      selection,
      targetCamera,
      target,
      fitDirection,
    );
  }
  targetCamera.updateMatrixWorld(true);
  return { camera: targetCamera, target };
}

function makeLineLayer(geometry, width, resolution, options) {
  const material = new LineMaterial({
    linewidth: width,
    vertexColors: options.vertexColors,
    color: options.color,
    transparent: options.transparent === undefined ? options.opacity < 1 : options.transparent,
    opacity: options.opacity,
    depthTest: options.depthTest === undefined ? true : options.depthTest,
    depthWrite: options.depthWrite,
    alphaToCoverage: true,
    clippingPlanes: options.clippingPlanes || null,
  });
  material.toneMapped = options.toneMapped === undefined ? true : options.toneMapped;
  material.resolution.copy(resolution);
  const line = new LineSegments2(geometry, material);
  line.renderOrder = options.renderOrder;
  return { line, material };
}

function polymerPaths(indices) {
  const paths = [];
  let current = null;
  for (let offset = 0; offset < indices.length; offset += 2) {
    const left = indices[offset];
    const right = indices[offset + 1];
    if (current && current[current.length - 1] === left) {
      current.push(right);
    } else {
      current = [left, right];
      paths.push(current);
    }
  }
  return paths;
}

function buildTubeGeometry(indices, colorSource, radius, radialSegments) {
  const paths = polymerPaths(indices);
  const ringCount = paths.reduce((total, path) => total + path.length, 0);
  const ringVertexCount = ringCount * radialSegments;
  const capVertexCount = paths.length * 2;
  const vertexCount = ringVertexCount + capVertexCount;
  const triangleIndexCount = (indices.length / 2) * radialSegments * 6
    + paths.length * radialSegments * 6;
  const vertexPositions = new Float32Array(vertexCount * 3);
  const vertexNormals = new Float32Array(vertexCount * 3);
  const vertexColors = new Float32Array(vertexCount * 3);
  const triangleIndices = new Uint32Array(triangleIndexCount);

  const center = new THREE.Vector3();
  const previous = new THREE.Vector3();
  const next = new THREE.Vector3();
  const tangent = new THREE.Vector3();
  const previousTangent = new THREE.Vector3();
  const normal = new THREE.Vector3();
  const binormal = new THREE.Vector3();
  const transportAxis = new THREE.Vector3();
  const helperAxis = new THREE.Vector3();
  const radialNormal = new THREE.Vector3();
  let ringOffset = 0;
  let indexOffset = 0;

  function readParticle(particleIndex, target) {
    return target.fromArray(positions, particleIndex * 3);
  }

  function writeVertex(vertexIndex, point, vertexNormal, particleIndex) {
    const target = vertexIndex * 3;
    vertexPositions[target] = point.x;
    vertexPositions[target + 1] = point.y;
    vertexPositions[target + 2] = point.z;
    vertexNormals[target] = vertexNormal.x;
    vertexNormals[target + 1] = vertexNormal.y;
    vertexNormals[target + 2] = vertexNormal.z;
    const colorOffset = particleIndex * 3;
    if (typeof colorSource === "function") {
      colorSource(particleIndex, vertexColors, target);
    } else {
      vertexColors[target] = colorSource[colorOffset];
      vertexColors[target + 1] = colorSource[colorOffset + 1];
      vertexColors[target + 2] = colorSource[colorOffset + 2];
    }
  }

  paths.forEach((path, pathIndex) => {
    const pathRingStart = ringOffset;
    for (let pointIndex = 0; pointIndex < path.length; pointIndex += 1) {
      const particleIndex = path[pointIndex];
      readParticle(particleIndex, center);
      if (pointIndex === 0) {
        readParticle(path[1], next);
        tangent.subVectors(next, center).normalize();
      } else if (pointIndex === path.length - 1) {
        readParticle(path[pointIndex - 1], previous);
        tangent.subVectors(center, previous).normalize();
      } else {
        readParticle(path[pointIndex - 1], previous);
        readParticle(path[pointIndex + 1], next);
        tangent.subVectors(next, previous).normalize();
      }

      if (pointIndex === 0) {
        if (Math.abs(tangent.z) < 0.85) helperAxis.set(0, 0, 1);
        else helperAxis.set(0, 1, 0);
        normal.crossVectors(helperAxis, tangent).normalize();
      } else {
        transportAxis.crossVectors(previousTangent, tangent);
        if (transportAxis.lengthSq() > 1e-10) {
          transportAxis.normalize();
          const angle = Math.acos(THREE.MathUtils.clamp(previousTangent.dot(tangent), -1, 1));
          normal.applyAxisAngle(transportAxis, angle);
        }
        normal.addScaledVector(tangent, -normal.dot(tangent)).normalize();
      }
      binormal.crossVectors(tangent, normal).normalize();
      normal.crossVectors(binormal, tangent).normalize();

      for (let segment = 0; segment < radialSegments; segment += 1) {
        const angle = (segment / radialSegments) * Math.PI * 2;
        radialNormal.copy(normal).multiplyScalar(Math.cos(angle));
        radialNormal.addScaledVector(binormal, Math.sin(angle)).normalize();
        const vertexIndex = (ringOffset + pointIndex) * radialSegments + segment;
        const point = next.copy(center).addScaledVector(radialNormal, radius);
        writeVertex(vertexIndex, point, radialNormal, particleIndex);
      }
      previousTangent.copy(tangent);
    }

    for (let pointIndex = 0; pointIndex < path.length - 1; pointIndex += 1) {
      const currentRing = (pathRingStart + pointIndex) * radialSegments;
      const nextRing = currentRing + radialSegments;
      for (let segment = 0; segment < radialSegments; segment += 1) {
        const following = (segment + 1) % radialSegments;
        const a = currentRing + segment;
        const b = currentRing + following;
        const c = nextRing + segment;
        const d = nextRing + following;
        triangleIndices[indexOffset++] = a;
        triangleIndices[indexOffset++] = b;
        triangleIndices[indexOffset++] = c;
        triangleIndices[indexOffset++] = b;
        triangleIndices[indexOffset++] = d;
        triangleIndices[indexOffset++] = c;
      }
    }

    const firstRing = pathRingStart * radialSegments;
    const lastRing = (pathRingStart + path.length - 1) * radialSegments;
    const startCap = ringVertexCount + pathIndex * 2;
    const endCap = startCap + 1;
    readParticle(path[0], center);
    readParticle(path[1], next);
    tangent.subVectors(next, center).normalize().negate();
    writeVertex(startCap, center, tangent, path[0]);
    readParticle(path[path.length - 1], center);
    readParticle(path[path.length - 2], previous);
    tangent.subVectors(center, previous).normalize();
    writeVertex(endCap, center, tangent, path[path.length - 1]);
    for (let segment = 0; segment < radialSegments; segment += 1) {
      const following = (segment + 1) % radialSegments;
      triangleIndices[indexOffset++] = startCap;
      triangleIndices[indexOffset++] = firstRing + following;
      triangleIndices[indexOffset++] = firstRing + segment;
      triangleIndices[indexOffset++] = endCap;
      triangleIndices[indexOffset++] = lastRing + segment;
      triangleIndices[indexOffset++] = lastRing + following;
    }
    ringOffset += path.length;
  });

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(vertexPositions, 3));
  geometry.setAttribute("normal", new THREE.BufferAttribute(vertexNormals, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(vertexColors, 3));
  geometry.setIndex(new THREE.BufferAttribute(triangleIndices, 1));
  geometry.computeBoundingBox();
  geometry.computeBoundingSphere();
  geometry.userData = {
    paths: paths.length,
    rings: ringCount,
    radialSegments,
    triangles: triangleIndices.length / 3,
  };
  return geometry;
}

function addStudioLights(targetScene, targetCamera, selection, profile) {
  const center = new THREE.Box3(
    new THREE.Vector3(...selection.bbox[0]),
    new THREE.Vector3(...selection.bbox[1]),
  ).getCenter(new THREE.Vector3());
  const extent = new THREE.Vector3(...selection.bbox[1]).sub(new THREE.Vector3(...selection.bbox[0]));
  const radius = Math.max(extent.length() * 0.5, 1);
  const forward = center.clone().sub(targetCamera.position).normalize();
  const right = new THREE.Vector3().crossVectors(forward, targetCamera.up).normalize();
  const up = new THREE.Vector3().crossVectors(right, forward).normalize();

  const hemisphere = new THREE.HemisphereLight(0xf8fbff, 0x263237, profile.lighting.hemisphere);
  const key = new THREE.DirectionalLight(0xfff9f0, profile.lighting.key);
  const fill = new THREE.DirectionalLight(0xc9ddff, profile.lighting.fill);
  const rim = new THREE.DirectionalLight(0xfff0dc, profile.lighting.rim);
  key.position.copy(center)
    .addScaledVector(right, -radius * 1.7)
    .addScaledVector(up, radius * 2.1)
    .addScaledVector(forward, -radius * 2.4);
  fill.position.copy(center)
    .addScaledVector(right, radius * 2.0)
    .addScaledVector(up, -radius * 0.8)
    .addScaledVector(forward, -radius * 1.3);
  rim.position.copy(center)
    .addScaledVector(right, radius * 0.6)
    .addScaledVector(up, radius * 1.4)
    .addScaledVector(forward, radius * 2.2);
  [key, fill, rim].forEach((light) => {
    light.target.position.copy(center);
    targetScene.add(light, light.target);
  });
  targetScene.add(hemisphere);
}

function makeTubeSilhouette(geometry, expansion, color, clippingPlanes = null, options = {}) {
  const material = new THREE.MeshBasicMaterial({
    color,
    side: THREE.BackSide,
    transparent: options.opacity !== undefined && options.opacity < 1,
    opacity: options.opacity === undefined ? 1 : options.opacity,
    depthTest: options.depthTest === undefined ? true : options.depthTest,
    depthWrite: options.depthWrite === undefined ? true : options.depthWrite,
    toneMapped: options.toneMapped === undefined ? true : options.toneMapped,
    clippingPlanes,
  });
  material.onBeforeCompile = (shader) => {
    shader.uniforms.outlineExpansion = { value: expansion };
    shader.vertexShader = shader.vertexShader
      .replace(
        "#include <common>",
        `#include <common>
        uniform float outlineExpansion;`,
      )
      .replace(
        "#include <begin_vertex>",
        "vec3 transformed = vec3( position ) + normal * outlineExpansion;",
      );
  };
  material.customProgramCacheKey = () => "tube-silhouette-v1";
  const mesh = new THREE.Mesh(geometry, material);
  mesh.renderOrder = 0;
  return { mesh, material };
}

function makeParticleMaterial(pointSize, clippingPlanes = null) {
  return new THREE.ShaderMaterial({
    uniforms: {
      pointSize: { value: pointSize },
      lightDirection: { value: new THREE.Vector3(-0.42, 0.56, 0.72).normalize() },
    },
    vertexShader: `
      #include <clipping_planes_pars_vertex>
      uniform float pointSize;
      attribute vec3 color;
      varying vec3 vColor;
      void main() {
        vColor = color;
        vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
        gl_Position = projectionMatrix * mvPosition;
        gl_PointSize = pointSize;
        #include <clipping_planes_vertex>
      }
    `,
    fragmentShader: `
      #include <clipping_planes_pars_fragment>
      uniform vec3 lightDirection;
      varying vec3 vColor;
      void main() {
        #include <clipping_planes_fragment>
        vec2 point = gl_PointCoord * 2.0 - 1.0;
        float radius2 = dot(point, point);
        if (radius2 > 1.0) discard;
        vec3 normal = normalize(vec3(point, sqrt(max(0.0, 1.0 - radius2))));
        float diffuse = 0.46 + 0.54 * max(dot(normal, lightDirection), 0.0);
        float specular = 0.20 * pow(max(dot(normal, normalize(vec3(-0.2, 0.35, 1.0))), 0.0), 18.0);
        float rim = smoothstep(0.72, 1.0, sqrt(radius2));
        vec3 shaded = vColor * diffuse + vec3(specular);
        shaded *= 1.0 - rim * 0.22;
        float alpha = 1.0 - smoothstep(0.90, 1.0, radius2);
        gl_FragColor = vec4(shaded, alpha);
      }
    `,
    transparent: true,
    depthTest: true,
    depthWrite: true,
    clipping: true,
    clippingPlanes,
  });
}

function sectionStatistics(state, selection) {
  if (!state.enabled) return null;
  const selected = selection.particleIndices;
  const particleCount = selected ? selected.length : payload.stats.particles;
  const signedDistance = (particle) => {
    const offset = particle * 3;
    return state.normal.x * (positions[offset] - state.center.x)
      + state.normal.y * (positions[offset + 1] - state.center.y)
      + state.normal.z * (positions[offset + 2] - state.center.z);
  };
  const half = state.thickness * 0.5;
  let visibleParticles = 0;
  for (let index = 0; index < particleCount; index += 1) {
    const particle = selected ? selected[index] : index;
    if (Math.abs(signedDistance(particle)) <= half) visibleParticles += 1;
  }
  let intersectingBonds = 0;
  for (let index = 0; index < selection.bondIndices.length; index += 2) {
    const left = signedDistance(selection.bondIndices[index]);
    const right = signedDistance(selection.bondIndices[index + 1]);
    if (Math.min(left, right) <= half && Math.max(left, right) >= -half) {
      intersectingBonds += 1;
    }
  }
  return {
    enabled: true,
    thickness: state.thickness,
    offset: state.offset,
    normal: state.normal.toArray(),
    center: state.center.toArray(),
    particles: visibleParticles,
    bonds: intersectingBonds,
  };
}

function annotationSectionBondCount(state, endpoints) {
  if (!state.enabled) return null;
  const half = state.thickness * 0.5;
  const signedDistance = (particle) => {
    const offset = particle * 3;
    return state.normal.x * (positions[offset] - state.center.x)
      + state.normal.y * (positions[offset + 1] - state.center.y)
      + state.normal.z * (positions[offset + 2] - state.center.z);
  };
  let count = 0;
  for (let index = 0; index < endpoints.length; index += 2) {
    const left = signedDistance(endpoints[index]);
    const right = signedDistance(endpoints[index + 1]);
    if (Math.min(left, right) <= half && Math.max(left, right) >= -half) count += 1;
  }
  return count;
}

function formatFeatureValue(value) {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return "";
  const absolute = Math.abs(numeric);
  if (absolute !== 0 && (absolute < 0.001 || absolute >= 10000)) return numeric.toExponential(2);
  return Number(numeric.toPrecision(4)).toString();
}

function compositeExportCanvas(exportCanvas, request) {
  const finalCanvas = document.createElement("canvas");
  finalCanvas.width = exportCanvas.width;
  finalCanvas.height = exportCanvas.height;
  const context = finalCanvas.getContext("2d");
  context.drawImage(exportCanvas, 0, 0);
  const mapping = payload.dataset.metadata && payload.dataset.metadata.color_mapping;
  if (!request.style.show_colorbar || !mapping || mapping.kind !== "continuous") {
    return { canvas: finalCanvas, colorbar: null };
  }

  const scale = Math.max(finalCanvas.width / 1400, 0.5);
  const barWidth = Math.round(190 * scale);
  const barHeight = Math.max(Math.round(14 * scale), 8);
  const margin = Math.round(42 * scale);
  const x = finalCanvas.width - margin - barWidth;
  const y = finalCanvas.height - margin - barHeight;
  const gradient = context.createLinearGradient(x, 0, x + barWidth, 0);
  mapping.stops.forEach(([position, color]) => gradient.addColorStop(position, color));
  context.fillStyle = gradient;
  context.fillRect(x, y, barWidth, barHeight);
  context.strokeStyle = "rgba(26, 35, 37, 0.72)";
  context.lineWidth = Math.max(1, scale);
  context.strokeRect(x, y, barWidth, barHeight);

  const fontSize = Math.max(Math.round(12 * scale), 10);
  context.font = `${fontSize}px Inter, Arial, sans-serif`;
  context.fillStyle = "#172126";
  context.textBaseline = "bottom";
  const baselineMode = request.style.profile === "section-pymol-v1";
  const minimum = baselineMode ? "Min" : formatFeatureValue(mapping.vmin);
  const maximum = baselineMode ? "Max" : formatFeatureValue(mapping.vmax);
  context.textAlign = "left";
  context.fillText(minimum, x, y - Math.round(5 * scale));
  context.textAlign = "right";
  context.fillText(maximum, x + barWidth, y - Math.round(5 * scale));
  if (!baselineMode) {
    context.textAlign = "left";
    context.fillStyle = "#34434a";
    context.fillText(mapping.label || mapping.field || "Feature", x, y - Math.round(24 * scale));
  }
  return {
    canvas: finalCanvas,
    colorbar: {
      label: mapping.label || mapping.field || null,
      cmap: mapping.cmap,
      vmin: mapping.vmin,
      vmax: mapping.vmax,
    },
  };
}

async function renderHighQuality(configOrProfile = null, options = {}) {
  const startedAt = performance.now();
  const request = normalizeRenderConfig(
    typeof configOrProfile === "string"
      ? legacyRenderConfig(configOrProfile, options)
      : (configOrProfile || interactiveRenderConfig()),
  );
  const profile = request.profile;
  const profileName = request.style.profile;
  const selected = selectionForRegions(request.view.regions);
  const selection = selected.selection;
  const sourceAspect = Math.max(renderer.domElement.clientWidth / renderer.domElement.clientHeight, 0.1);
  const longEdge = request.raster.long_edge;
  const width = request.raster.width || (sourceAspect >= 1 ? longEdge : Math.round(longEdge * sourceAspect));
  const height = request.raster.height || (sourceAspect >= 1 ? Math.round(longEdge / sourceAspect) : longEdge);
  const exportCanvas = document.createElement("canvas");
  const exportRenderer = new THREE.WebGLRenderer({
    canvas: exportCanvas,
    antialias: true,
    alpha: false,
    powerPreference: "default",
    preserveDrawingBuffer: true,
  });
  exportRenderer.setPixelRatio(1);
  exportRenderer.setSize(width, height, false);
  exportRenderer.setClearColor(request.style.background, 1);
  exportRenderer.outputColorSpace = THREE.SRGBColorSpace;
  exportRenderer.toneMapping = THREE.ACESFilmicToneMapping;
  exportRenderer.toneMappingExposure = profile.exposure;
  await new Promise((resolve) => requestAnimationFrame(resolve));

  const exportScene = new THREE.Scene();
  const configured = configuredCamera(request.view.camera, selection, width / height);
  const exportCamera = configured.camera;
  const exportTarget = configured.target;
  const exportSection = sectionState(request.view.section, exportCamera, exportTarget);
  const clippingPlanes = exportSection.enabled ? exportSection.planes : null;
  exportRenderer.localClippingEnabled = exportSection.enabled;
  const resolution = new THREE.Vector2(width, height);
  const profileScale = width / 1400;
  const exportColors = profileColorArray(profileName, exportCamera, selection);
  const exportGeometries = [];
  const exportMaterials = [];
  let renderPasses = 0;
  let geometryStats = null;
  let studioLightsAdded = false;
  const annotationReport = [];

  const addParticles = () => {
    const exportParticleGeometry = new THREE.BufferGeometry();
    exportParticleGeometry.setAttribute("position", positionAttribute);
    exportParticleGeometry.setAttribute("color", new THREE.BufferAttribute(exportColors, 3));
    if (selection.particleIndices) {
      exportParticleGeometry.setIndex(new THREE.BufferAttribute(selection.particleIndices, 1));
    }
    const exportParticleMaterial = makeParticleMaterial(
      request.style.particle_size * profileScale,
      clippingPlanes,
    );
    const exportParticles = new THREE.Points(exportParticleGeometry, exportParticleMaterial);
    exportParticles.renderOrder = 3;
    exportScene.add(exportParticles);
    exportGeometries.push(exportParticleGeometry);
    exportMaterials.push(exportParticleMaterial);
    renderPasses += 1;
  };

  const hasBonds = request.layers.show_bonds && selection.bondIndices.length > 0;
  const hasRequestedAnnotations = request.layers.show_annotations
    && annotationLayers.some((layer) => layer.visible && layer.annotations.length > 0);
  if (!hasBonds && !hasRequestedAnnotations && !request.particlesSpecified) {
    request.layers.show_particles = true;
  }
  if (hasBonds && profile.mode.startsWith("layered")) {
    const arrays = lineArrays(selection.bondIndices, exportColors, false);
    const exportLineGeometry = new LineSegmentsGeometry();
    exportLineGeometry.setPositions(arrays.linePositions);
    exportLineGeometry.setColors(arrays.lineColors);
    exportGeometries.push(exportLineGeometry);
    const outline = request.style.silhouette_width > request.style.tube_width
      ? makeLineLayer(exportLineGeometry, request.style.silhouette_width * profileScale, resolution, {
        vertexColors: false,
        color: profile.outline,
        opacity: profile.outlineOpacity,
        depthWrite: true,
        renderOrder: 0,
        clippingPlanes,
      })
      : null;
    const inner = makeLineLayer(exportLineGeometry, request.style.tube_width * profileScale, resolution, {
      vertexColors: true,
      color: 0xffffff,
      opacity: 0.96,
      depthWrite: false,
      renderOrder: 1,
      clippingPlanes,
    });
    const highlight = profile.highlightOpacity > 0
      ? makeLineLayer(exportLineGeometry, profile.highlightWidth * profileScale, resolution, {
        vertexColors: false,
        color: 0xffffff,
        opacity: profile.highlightOpacity,
        depthWrite: false,
        renderOrder: 2,
        clippingPlanes,
      })
      : null;
    if (outline) {
      exportScene.add(outline.line);
      exportMaterials.push(outline.material);
      renderPasses += 1;
    }
    exportScene.add(inner.line);
    exportMaterials.push(inner.material);
    renderPasses += 1;
    if (highlight) {
      exportScene.add(highlight.line);
      exportMaterials.push(highlight.material);
      renderPasses += 1;
    }
  } else if (hasBonds) {
    const targetDistance = exportCamera.position.distanceTo(exportTarget);
    const worldHeight = 2 * targetDistance * Math.tan(THREE.MathUtils.degToRad(exportCamera.fov) / 2);
    const worldPerPixel = worldHeight / height;
    const tubeDiameter = request.style.tube_width * profileScale * worldPerPixel;
    const tubeGeometry = buildTubeGeometry(
      selection.bondIndices,
      exportColors,
      tubeDiameter * 0.5,
      profile.radialSegments,
    );
    exportGeometries.push(tubeGeometry);
    geometryStats = tubeGeometry.userData;

    const outlineExpansion = Math.max(
      0,
      (request.style.silhouette_width - request.style.tube_width) * profileScale * worldPerPixel * 0.5,
    );
    if (outlineExpansion > 0) {
      const silhouette = makeTubeSilhouette(
        tubeGeometry,
        outlineExpansion,
        profile.silhouette,
        clippingPlanes,
      );
      exportScene.add(silhouette.mesh);
      exportMaterials.push(silhouette.material);
      renderPasses += 1;
    }

    const tubeMaterial = new THREE.MeshStandardMaterial({
      vertexColors: true,
      roughness: profile.lighting.roughness,
      metalness: 0,
      side: THREE.FrontSide,
      clippingPlanes,
    });
    tubeMaterial.dithering = true;
    const tubeMesh = new THREE.Mesh(tubeGeometry, tubeMaterial);
    tubeMesh.renderOrder = 1;
    exportScene.add(tubeMesh);
    exportMaterials.push(tubeMaterial);
    addStudioLights(exportScene, exportCamera, selection, profile);
    studioLightsAdded = true;
    renderPasses += 1;
  }

  if (request.layers.show_particles) addParticles();

  annotationLayers.forEach((layer, layerIndex) => {
    const layerVisible = request.layers.show_annotations && layer.visible;
    const layerReport = {
      id: layer.id,
      label: layer.label,
      visible: layerVisible,
      order: layer.order,
      annotations: [],
    };
    layer.annotations.forEach((annotation, annotationIndex) => {
      const endpoints = layerVisible
        ? annotationBondEndpoints(annotation, selection)
        : new Uint32Array();
      const renderedParticles = new Set(endpoints).size;
      const report = {
        id: annotation.id,
        label: annotation.label,
        order: annotation.order,
        target: {
          kind: annotation.target.kind,
          selector: annotation.target.selector,
          chromosomes: annotation.target.chromosomes,
          resolved_particles: annotation.target.particle_count,
          resolved_bonds: annotation.target.bond_count,
          rendered_particles: renderedParticles,
          rendered_bonds: endpoints.length / 2,
          section_bonds: annotationSectionBondCount(exportSection, endpoints),
          anchor: annotation.target.anchor,
          bbox: annotation.target.bbox,
        },
        marks: [],
      };
      annotation.marks.forEach((mark, markIndex) => {
        if (mark.kind !== "tube") throw new Error(`Unsupported annotation mark: ${mark.kind}`);
        report.marks.push({
          kind: mark.kind,
          requested_color: mark.requested_color,
          resolved_color: mark.color,
          tube_width: mark.tube_width,
          outline_width: mark.outline_width,
          requested_outline_color: mark.requested_outline_color,
          resolved_outline_color: mark.outline_color,
          material: mark.material,
          opacity: mark.opacity,
          occlusion: mark.occlusion,
          profile_color_transforms_bypassed: mark.material === "unlit",
          tone_mapping_bypassed: mark.material === "unlit",
          outline_tone_mapping_bypassed: true,
        });
        if (!endpoints.length) return;

        const targetDistance = exportCamera.position.distanceTo(exportTarget);
        const worldHeight = 2 * targetDistance * Math.tan(THREE.MathUtils.degToRad(exportCamera.fov) / 2);
        const worldPerPixel = worldHeight / height;
        const radius = mark.tube_width * profileScale * worldPerPixel * 0.5;
        const colorWriter = annotationColorWriter(mark, profile, exportCamera, selection);
        const geometry = buildTubeGeometry(endpoints, colorWriter, radius, profile.radialSegments || 8);
        exportGeometries.push(geometry);
        const depthTest = mark.occlusion === "depth_tested";
        const depthWrite = depthTest && mark.opacity === 1;
        const renderOrder = mark.occlusion === "always_visible"
          ? 10_000 + layerIndex * 100 + annotationIndex * 4 + markIndex * 2
          : 100 + layerIndex * 100 + annotationIndex * 4 + markIndex * 2;
        const outlineExpansion = Math.max(
          0,
          (mark.outline_width - mark.tube_width) * profileScale * worldPerPixel * 0.5,
        );
        if (outlineExpansion > 0 && mark.outline_color) {
          const outline = makeTubeSilhouette(
            geometry,
            outlineExpansion,
            mark.outline_color,
            clippingPlanes,
            {
              opacity: mark.opacity,
              depthTest,
              depthWrite,
              toneMapped: false,
            },
          );
          outline.mesh.renderOrder = renderOrder;
          exportScene.add(outline.mesh);
          exportMaterials.push(outline.material);
          renderPasses += 1;
        }

        let material;
        if (mark.material === "unlit" || !profile.lighting) {
          material = new THREE.MeshBasicMaterial({
            vertexColors: true,
            side: THREE.FrontSide,
            transparent: mark.opacity < 1,
            opacity: mark.opacity,
            depthTest,
            depthWrite,
            toneMapped: mark.material !== "unlit",
            clippingPlanes,
          });
        } else {
          material = new THREE.MeshStandardMaterial({
            vertexColors: true,
            roughness: profile.lighting.roughness,
            metalness: 0,
            side: THREE.FrontSide,
            transparent: mark.opacity < 1,
            opacity: mark.opacity,
            depthTest,
            depthWrite,
            clippingPlanes,
          });
          if (!studioLightsAdded) {
            addStudioLights(exportScene, exportCamera, selection, profile);
            studioLightsAdded = true;
          }
        }
        const mesh = new THREE.Mesh(geometry, material);
        mesh.renderOrder = renderOrder + 1;
        exportScene.add(mesh);
        exportMaterials.push(material);
        renderPasses += 1;
      });
      layerReport.annotations.push(report);
    });
    annotationReport.push(layerReport);
  });

  if (renderPasses === 0) {
    throw new Error("Render request has no visible particles, polymer bonds, or annotations");
  }

  await new Promise((resolve) => requestAnimationFrame(resolve));
  exportRenderer.render(exportScene, exportCamera);
  exportRenderer.render(exportScene, exportCamera);
  const composited = compositeExportCanvas(exportCanvas, request);
  const blob = await new Promise((resolve, reject) => {
    composited.canvas.toBlob((result) => {
      if (result) resolve(result);
      else reject(new Error("PNG encoder returned no data"));
    }, "image/png");
  });

  exportGeometries.forEach((geometry) => geometry.dispose());
  exportMaterials.forEach((material) => material.dispose());
  exportRenderer.dispose();
  exportRenderer.forceContextLoss();
  return {
    blob,
    width,
    height,
    profile: profile.label,
    renderPasses,
    geometryStats,
    colorbar: composited.colorbar,
    request: {
      schema_version: request.schema_version,
      view: request.view,
      layers: request.layers,
      style: request.style,
      raster: request.raster,
    },
    selection: {
      regions: selected.regions.map((region) => region.raw),
      particles: selection.particleCount,
      bonds: selection.bondCount,
    },
    annotations: annotationReport,
    section: sectionStatistics(exportSection, selection),
    camera: {
      mode: request.view.camera.mode,
      position: exportCamera.position.toArray(),
      target: exportTarget.toArray(),
      up: exportCamera.up.toArray(),
      fov: exportCamera.fov,
    },
    sceneSha256: payload.dataset.sha256,
    renderer: {
      threeRevision: THREE.REVISION,
      software: softwareRenderer,
      label: rendererLabel,
    },
    renderTimeMs: Math.round(performance.now() - startedAt),
  };
}

function triggerDownload(blob, filename) {
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = filename;
  link.click();
  setTimeout(() => URL.revokeObjectURL(link.href), 1000);
}

async function downloadPng(config = null) {
  saveButton.disabled = true;
  const previousLabel = saveButton.textContent;
  saveButton.textContent = "Rendering";
  fpsStatus.textContent = "Rendering PNG";
  try {
    const result = await renderHighQuality(config);
    const stem = payload.dataset.name.replace(/\.[^.]+$/, "");
    triggerDownload(result.blob, `${stem}.${result.profile}.${result.width}x${result.height}.png`);
    return {
      width: result.width,
      height: result.height,
      profile: result.profile,
      renderPasses: result.renderPasses,
      geometryStats: result.geometryStats,
      size: result.blob.size,
      request: result.request,
      selection: result.selection,
      annotations: result.annotations,
      section: result.section,
      colorbar: result.colorbar,
      camera: result.camera,
      sceneSha256: result.sceneSha256,
      renderer: result.renderer,
      renderTimeMs: result.renderTimeMs,
    };
  } finally {
    saveButton.disabled = false;
    saveButton.textContent = previousLabel;
  }
}

saveButton.addEventListener("click", () => {
  downloadPng().catch((error) => {
    errorState.style.display = "grid";
    errorState.textContent = `PNG export failed: ${error.message}`;
  });
});

controls.addEventListener("start", () => {
  controlsActive = true;
  tooltip.style.display = "none";
});
controls.addEventListener("end", () => {
  controlsActive = false;
});
controls.addEventListener("change", () => {
  if (sectionToggle.checked) updateMonitorSectionPlanes();
});

const raycaster = new THREE.Raycaster();
raycaster.params.Points.threshold = payload.viewer.base_point_size * 1.7;
const pointer = new THREE.Vector2();

function particleLabel(index) {
  const chrom = payload.chromosomes[chromosomeCodes[index]];
  const position = genomicPositions[index].toLocaleString();
  return { chrom, position, text: `${chrom}:${position}` };
}

function pickParticle(event, force = false) {
  if (controlsActive || !particles.visible) return null;
  const now = performance.now();
  if (!force && now - lastPickTime < 90) return null;
  lastPickTime = now;

  const rect = renderer.domElement.getBoundingClientRect();
  pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
  raycaster.setFromCamera(pointer, camera);
  const hit = raycaster.intersectObject(particles, false)[0];
  if (!hit) {
    tooltip.style.display = "none";
    return null;
  }

  const label = particleLabel(hit.index);
  tooltip.replaceChildren();
  const title = document.createElement("strong");
  title.textContent = label.chrom;
  const detail = document.createElement("span");
  detail.textContent = `${label.position} bp`;
  tooltip.append(title, detail);
  tooltip.style.left = `${Math.min(event.clientX + 14, window.innerWidth - 190)}px`;
  tooltip.style.top = `${Math.min(event.clientY + 14, window.innerHeight - 76)}px`;
  tooltip.style.display = "block";
  return { ...hit, label };
}

renderer.domElement.addEventListener("pointermove", (event) => pickParticle(event));
renderer.domElement.addEventListener("pointerleave", () => {
  tooltip.style.display = "none";
});
renderer.domElement.addEventListener("click", (event) => {
  const hit = pickParticle(event, true);
  if (!hit) return;
  selectionHalo.position.fromArray(positions, hit.index * 3);
  selectionHalo.visible = true;
  selectionStatus.textContent = `Selected ${hit.label.text} | particle ${hit.index.toLocaleString()}`;
});

function positionRegionPanel() {
  const topbar = document.querySelector(".topbar");
  regionPanel.style.top = `${Math.ceil(topbar.getBoundingClientRect().bottom + 8)}px`;
}

function resize() {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
  if (wideBondMaterial) wideBondMaterial.resolution.set(window.innerWidth, window.innerHeight);
  if (hasAnnotations) rebuildMonitorAnnotations();
  if (!regionPanel.hidden) positionRegionPanel();
}
window.addEventListener("resize", resize);

function animate(now) {
  const delta = now - lastFrameTime;
  lastFrameTime = now;
  controls.update(delta / 1000);
  renderer.render(scene, camera);

  fpsFrames += 1;
  if (now - fpsWindowStart >= 500) {
    const fps = (fpsFrames * 1000) / (now - fpsWindowStart);
    fpsStatus.textContent = `${fps.toFixed(0)} FPS`;
    fpsFrames = 0;
    fpsWindowStart = now;
    const rendererMode = softwareRenderer ? "Software WebGL 2" : "WebGL 2";
    rendererStatus.textContent = `Three.js r${THREE.REVISION} | ${rendererMode} | ${renderer.info.render.calls} draw calls`;
  }

  if (firstFrameTime === null) {
    firstFrameTime = now;
    document.body.dataset.ready = "1";
    document.body.dataset.firstFrameMs = String(Math.round(now - payload.generated_at_ms));
  }
}

renderRegionList();
applySelection(allSelection, false);
if (payload.viewer.initial_regions && payload.viewer.initial_regions.length) {
  try {
    addRegions(payload.viewer.initial_regions.join(";"));
  } catch (error) {
    regionError.textContent = error.message;
    applySelection(allSelection, false);
  }
}
particles.visible = particlesToggle.checked;
bonds.visible = bondsToggle.checked;
monitorAnnotations.visible = annotationsToggle.checked;
const initialCameraConfig = embeddedRenderConfig.view.camera || { mode: "fit", fov: 34, up: [0, 1, 0] };
if (initialCameraConfig.mode === "explicit") {
  const configured = configuredCamera(initialCameraConfig, activeSelection, camera.aspect);
  camera.copy(configured.camera);
  controls.target.copy(configured.target);
  controls.update();
} else {
  camera.fov = initialCameraConfig.fov || 34;
  camera.up.fromArray(initialCameraConfig.up || [0, 1, 0]);
  if (initialCameraConfig.direction) {
    frameCamera(
      activeSelection,
      camera,
      controls.target,
      new THREE.Vector3(...initialCameraConfig.direction).normalize(),
    );
    controls.update();
  } else {
    fitView(activeSelection);
  }
}
configureSectionInputs(activeSelection);
updateSectionControls();
rendererStatus.textContent = `Three.js r${THREE.REVISION} | ${softwareRenderer ? "Software WebGL 2" : "WebGL 2"}`;
renderer.setAnimationLoop(animate);

window.__demo = {
  scene,
  camera,
  renderer,
  controls,
  particles,
  bonds,
  payload,
  fitView,
  parseRegion,
  addRegions,
  renderHighQuality,
  downloadPng,
  clearRegions() {
    activeRegions = [];
    applyRegions();
  },
  getState() {
    return {
      ready: document.body.dataset.ready === "1",
      camera: camera.position.toArray(),
      target: controls.target.toArray(),
      particleVisible: particles.visible,
      bondVisible: bonds.visible,
      annotationVisible: monitorAnnotations.visible,
      pointDrawCount: activeSelection.particleCount,
      bondDrawCount: activeSelection.bondCount,
      regionCount: activeRegions.length,
      regionLabels: activeRegions.map((region) => region.raw),
      bondSize: Number(bondSizeInput.value),
      section: {
        enabled: sectionToggle.checked,
        thickness: Number(sectionThicknessInput.value),
        offset: Number(sectionOffsetInput.value),
      },
      wideBondMode: Boolean(wideBonds && wideBonds.visible),
      background: `#${renderer.getClearColor(new THREE.Color()).getHexString()}`,
      drawCalls: renderer.info.render.calls,
      triangles: renderer.info.render.triangles,
      points: renderer.info.render.points,
      lines: renderer.info.render.lines,
      softwareRenderer,
      rendererLabel,
    };
  },
};

window.__hicBasicRenderer = {
  protocolVersion: RENDERER_PROTOCOL_VERSION,
  renderConfigSchemaVersion: 1,
  renderHighQuality,
  downloadPng,
  getState: window.__demo.getState,
};
