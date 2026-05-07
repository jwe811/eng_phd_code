import React from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { api, UofsObject } from '../lib/api';

const colors = [0x0B6A41, 0xDBCC52, 0x000000, 0x7a1f1f];
const colorLabels = ['#0B6A41', '#DBCC52', '#000000', '#7A1F1F'];

function labelSprite(text: string, color = '#000000') {
  const canvas = document.createElement('canvas');
  canvas.width = 256;
  canvas.height = 128;
  const context = canvas.getContext('2d')!;
  context.clearRect(0, 0, canvas.width, canvas.height);
  context.fillStyle = color;
  context.font = 'bold 80px sans-serif';
  context.textAlign = 'center';
  context.textBaseline = 'middle';
  context.fillText(text, canvas.width / 2, canvas.height / 2);
  const texture = new THREE.CanvasTexture(canvas);
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, transparent: true }));
  sprite.scale.set(1.4, 0.7, 1);
  return sprite;
}

function edgeCylinder(start: THREE.Vector3, end: THREE.Vector3, color: number, radius: number) {
  const midpoint = start.clone().add(end).multiplyScalar(0.5);
  const direction = end.clone().sub(start);
  const length = direction.length();
  const geometry = new THREE.CylinderGeometry(radius, radius, length, 12);
  const material = new THREE.MeshStandardMaterial({ color });
  const mesh = new THREE.Mesh(geometry, material);
  mesh.position.copy(midpoint);
  mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.normalize());
  return mesh;
}

export function PolygonPreview(props: { path: string; index: number }) {
  const mountRef = React.useRef<HTMLDivElement | null>(null);
  const controlsRef = React.useRef<OrbitControls | null>(null);
  const [data, setData] = React.useState<UofsObject | null>(null);
  const [error, setError] = React.useState('');

  React.useEffect(() => {
    let cancelled = false;
    setError('');
    setData(null);
    api.uofsObject(props.path, props.index)
      .then((object) => {
        if (!cancelled) setData(object);
      })
      .catch((err) => {
        if (!cancelled) setError(err instanceof Error ? err.message : String(err));
      });
    return () => {
      cancelled = true;
    };
  }, [props.path, props.index]);

  React.useEffect(() => {
    const mount = mountRef.current;
    if (!mount || !data) return;
    mount.innerHTML = '';
    const width = mount.clientWidth || 720;
    const height = mount.clientHeight || 420;
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0xf8fafc);
    const camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 1000);
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(width, height);
    mount.appendChild(renderer.domElement);

    const allPoints = data.polygons.flatMap((poly) => poly.points);
    const vectors = allPoints.map((point) => new THREE.Vector3(point[0], point[1], point[2]));
    const center = new THREE.Vector3();
    let bounds: THREE.Box3;

    if (data.L !== null && data.M !== null && data.span !== null) {
      bounds = new THREE.Box3(
        new THREE.Vector3(0, 0, 0),
        new THREE.Vector3(data.span, data.L, data.M)
      );
      center.set(data.span / 2, data.L / 2, data.M / 2);
    } else {
      for (const point of vectors) center.add(point);
      center.divideScalar(Math.max(1, vectors.length));
      bounds = new THREE.Box3().setFromPoints(vectors).expandByScalar(0.25);
    }

    const maxRadius = Math.max(
      4,
      ...vectors.map((point) => point.distanceTo(center)),
      ...(data.L !== null && data.M !== null && data.span !== null ? [data.L / 2, data.M / 2, data.span / 2] : [])
    );

    const initialPosition = new THREE.Vector3(center.x + maxRadius * 2, center.y + maxRadius * 2, center.z + maxRadius * 1.5);
    camera.position.copy(initialPosition);
    camera.up.set(0, 0, 1);
    camera.lookAt(center);

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.target.copy(center);
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controlsRef.current = controls;

    scene.add(new THREE.AmbientLight(0xffffff, 0.7));
    const light = new THREE.DirectionalLight(0xffffff, 0.7);
    light.position.set(4, 5, 7);
    scene.add(light);

    const gridDivs = Math.max(10, Math.ceil(maxRadius * 4));
    const finalDivs = gridDivs % 2 === 0 ? gridDivs : gridDivs + 1;
    const grid = new THREE.GridHelper(finalDivs, finalDivs, 0xd1d5db, 0xe5e7eb);
    grid.rotation.x = Math.PI / 2;
    grid.position.set(Math.round(center.x), Math.round(center.y), bounds.min.z);
    scene.add(grid);

    const size = bounds.getSize(new THREE.Vector3());
    const lengthX = size.x + 2;
    const lengthY = size.y + 2;
    const lengthZ = size.z + 2;

    scene.add(new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), bounds.min, lengthX, 0x000000, 0.35, 0.22));
    scene.add(new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), bounds.min, lengthY, 0x0B6A41, 0.35, 0.22));
    scene.add(new THREE.ArrowHelper(new THREE.Vector3(0, 0, 1), bounds.min, lengthZ, 0xDBCC52, 0.35, 0.22));
    const xLabel = labelSprite('x');
    xLabel.position.set(bounds.min.x + lengthX + 0.6, bounds.min.y, bounds.min.z);
    const yLabel = labelSprite('y', '#0B6A41');
    yLabel.position.set(bounds.min.x, bounds.min.y + lengthY + 0.6, bounds.min.z);
    const zLabel = labelSprite('z', '#000000');
    zLabel.position.set(bounds.min.x, bounds.min.y, bounds.min.z + lengthZ + 0.6);
    scene.add(xLabel, yLabel, zLabel);

    data.polygons.forEach((poly, index) => {
      const points = poly.points.map((point) => new THREE.Vector3(point[0], point[1], point[2]));
      const color = colors[index % colors.length];
      for (let pointIndex = 0; pointIndex < points.length - 1; pointIndex += 1) {
        scene.add(edgeCylinder(points[pointIndex], points[pointIndex + 1], color, 0.018 * maxRadius));
      }
      const nodeMaterial = new THREE.MeshStandardMaterial({ color: colors[index % colors.length] });
      for (const point of points) {
        const node = new THREE.Mesh(new THREE.SphereGeometry(0.0495 * maxRadius, 16, 12), nodeMaterial);
        node.position.copy(point);
        scene.add(node);
      }
    });

    let animationId: number;
    function animate() {
      animationId = requestAnimationFrame(animate);
      controls.update();
      renderer.render(scene, camera);
    }
    animate();

    return () => {
      cancelAnimationFrame(animationId);
      controls.dispose();
      renderer.dispose();
      mount.innerHTML = '';
    };
  }, [data]);

  return (
    <div className="visualizer">
      {error && <div className="errors">{error}</div>}
      {data?.validation_errors.length ? <div className="errors">{data.validation_errors.join(' ')}</div> : null}
      {!data && !error && <div className="canvas-wrap loading-canvas"><span>Loading object geometry...</span></div>}
      {data && (
        <>
          <div className="visualizer-toolbar">
            <div className="metadata-grid compact">
              <div><span>Object</span><strong>{data.index + 1} / {data.object_count}</strong></div>
              <div><span>Span (x)</span><strong>{data.span_x}</strong></div>
              <div><span>Contacts</span><strong>{data.contacts}</strong></div>
              <div><span>Shrink label</span><strong>{data.shrink_labels.join(', ')}</strong></div>
              {data.linking_number !== null && (
                <div><span>Linking #</span><strong>{data.linking_number}</strong></div>
              )}
              <div><span>Edges</span><strong>{data.polygons.reduce((sum, poly) => sum + poly.length, 0)}</strong></div>
            </div>
            <div className="button-strip">
              <button onClick={() => {
                if (controlsRef.current) {
                  controlsRef.current.reset();
                }
              }}>Reset view</button>
            </div>
          </div>
          <div className="legend">
            {data.polygons.map((poly, index) => (
              <span key={index}><i style={{ background: colorLabels[index % colorLabels.length] }} /> Polygon {index + 1} · {poly.length} edges</span>
            ))}
          </div>
          <div ref={mountRef} className="canvas-wrap" aria-label="3D polygon visualization" />
        </>
      )}
    </div>
  );
}
