import React from 'react';
import * as THREE from 'three';
import { api, UofsObject } from '../lib/api';

const colors = [0x0B6A41, 0xDBCC52, 0x000000, 0x7a1f1f];
const colorLabels = ['#0B6A41', '#DBCC52', '#000000', '#7A1F1F'];

function labelSprite(text: string, color = '#000000') {
  const canvas = document.createElement('canvas');
  canvas.width = 128;
  canvas.height = 64;
  const context = canvas.getContext('2d')!;
  context.fillStyle = 'rgba(255,255,255,0.92)';
  context.fillRect(0, 0, canvas.width, canvas.height);
  context.fillStyle = color;
  context.font = 'bold 34px sans-serif';
  context.textAlign = 'center';
  context.textBaseline = 'middle';
  context.fillText(text, canvas.width / 2, canvas.height / 2);
  const texture = new THREE.CanvasTexture(canvas);
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture }));
  sprite.scale.set(0.7, 0.35, 1);
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
  const resetRef = React.useRef<() => void>(() => undefined);
  const zoomRef = React.useRef<(factor: number) => void>(() => undefined);
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
    for (const point of vectors) center.add(point);
    center.divideScalar(Math.max(1, vectors.length));
    const maxRadius = Math.max(
      4,
      ...vectors.map((point) => point.distanceTo(center))
    );

    const initialPosition = new THREE.Vector3(center.x + maxRadius * 2.2, center.y + maxRadius * 2.2, center.z + maxRadius * 1.7);
    function resetCamera() {
      camera.position.copy(initialPosition);
      camera.lookAt(center);
      renderer.render(scene, camera);
    }
    function zoom(factor: number) {
      const offset = camera.position.clone().sub(center).multiplyScalar(factor);
      camera.position.copy(center.clone().add(offset));
      camera.lookAt(center);
      renderer.render(scene, camera);
    }
    resetCamera();
    resetRef.current = resetCamera;
    zoomRef.current = zoom;
    scene.add(new THREE.AmbientLight(0xffffff, 0.7));
    const light = new THREE.DirectionalLight(0xffffff, 0.7);
    light.position.set(4, 5, 7);
    scene.add(light);

    const grid = new THREE.GridHelper(Math.max(4, maxRadius * 3), Math.max(4, Math.ceil(maxRadius * 3)));
    grid.position.set(center.x, center.y, center.z - 0.02);
    scene.add(grid);

    const bounds = new THREE.Box3().setFromPoints(vectors).expandByScalar(0.25);
    scene.add(new THREE.Box3Helper(bounds, 0x0B6A41));

    const size = bounds.getSize(new THREE.Vector3());
    const axisLength = Math.max(1.5, Math.max(size.x, size.y, size.z) + 1);
    scene.add(new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), bounds.min, axisLength, 0x000000));
    scene.add(new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), bounds.min, axisLength, 0x0B6A41));
    scene.add(new THREE.ArrowHelper(new THREE.Vector3(0, 0, 1), bounds.min, axisLength, 0xDBCC52));
    const xLabel = labelSprite('x');
    xLabel.position.set(bounds.min.x + axisLength + 0.3, bounds.min.y, bounds.min.z);
    const yLabel = labelSprite('y', '#0B6A41');
    yLabel.position.set(bounds.min.x, bounds.min.y + axisLength + 0.3, bounds.min.z);
    const zLabel = labelSprite('z', '#000000');
    zLabel.position.set(bounds.min.x, bounds.min.y, bounds.min.z + axisLength + 0.3);
    scene.add(xLabel, yLabel, zLabel);

    data.polygons.forEach((poly, index) => {
      const points = poly.points.map((point) => new THREE.Vector3(point[0], point[1], point[2]));
      const color = colors[index % colors.length];
      for (let pointIndex = 0; pointIndex < points.length - 1; pointIndex += 1) {
        scene.add(edgeCylinder(points[pointIndex], points[pointIndex + 1], color, 0.02 * maxRadius));
      }
      const nodeMaterial = new THREE.MeshStandardMaterial({ color: colors[index % colors.length] });
      for (const point of points) {
        const node = new THREE.Mesh(new THREE.SphereGeometry(0.055 * maxRadius, 16, 12), nodeMaterial);
        node.position.copy(point);
        scene.add(node);
      }
    });

    let dragging = false;
    let lastX = 0;
    function rotate(delta: number) {
      const offset = camera.position.clone().sub(center);
      offset.applyAxisAngle(new THREE.Vector3(0, 0, 1), delta);
      camera.position.copy(center.clone().add(offset));
      camera.lookAt(center);
    }
    function onPointerDown(event: PointerEvent) {
      dragging = true;
      lastX = event.clientX;
    }
    function onPointerMove(event: PointerEvent) {
      if (!dragging) return;
      rotate((event.clientX - lastX) * 0.01);
      lastX = event.clientX;
      renderer.render(scene, camera);
    }
    function onPointerUp() {
      dragging = false;
    }
    renderer.domElement.addEventListener('pointerdown', onPointerDown);
    renderer.domElement.addEventListener('pointermove', onPointerMove);
    window.addEventListener('pointerup', onPointerUp);
    renderer.render(scene, camera);

    return () => {
      renderer.domElement.removeEventListener('pointerdown', onPointerDown);
      renderer.domElement.removeEventListener('pointermove', onPointerMove);
      window.removeEventListener('pointerup', onPointerUp);
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
              <div><span>Object</span><strong>{data.index + 1} of {data.object_count}</strong></div>
              <div><span>Polygons</span><strong>{data.polygons.length}</strong></div>
              <div><span>Total edges</span><strong>{data.polygons.reduce((sum, poly) => sum + poly.length, 0)}</strong></div>
            </div>
            <div className="button-strip">
              <button onClick={() => resetRef.current()}>Reset camera</button>
              <button onClick={() => zoomRef.current(0.8)}>Zoom in</button>
              <button onClick={() => zoomRef.current(1.25)}>Zoom out</button>
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
